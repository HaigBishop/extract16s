"""
asvs2truncspec_03_agg.py - Process HMMER outputs and generate .truncspec file

This script is Step 3 of the asvs2truncspec pipeline. It takes the HMMER outputs 
from Step 2 and converts them into a .truncspec file that specifies regions to 
extract from 16S sequences.

See asvs2truncspec-plan.md for the overall process and pipeline details.

Processing steps:
1. Parse Stockholm alignments to create model→reference mappings
2. Parse tblout files and apply hit filtering (E-value, coverage)
3. Select best hit per ASV per HMM, then apply best-HMM gate
4. Convert HMM coordinates to reference coordinates
5. Compute consensus regions per dataset using median + outlier filtering
6. Optional cross-domain bootstrapping to fill missing coords
7. Optional region redundancy minimisation to merge similar regions
8. Optional region renaming to simplify region names
9. Write outputs:
   - Final .truncspec file
   - Model→reference mapping TSVs
   - Per-ASV hit TSVs
   - Dataset region calls and QC summaries
   - Filter statistics
   - Region renaming map (if enabled)

Run from the extract16s directory:
    python Scripts/asvs2truncspec_03_agg.py
"""


# Imports
import os
import gzip
import datetime
import copy
import matplotlib.pyplot as plt


# =============================================================================
# CONFIGURATION - Modify these values as needed
# =============================================================================

# Output directory for info and intermediates (same as Step 1)
INFO_OUT_DIR = "/home/haig/Repos/micro16s/extract16s/asvs2truncspec_out/"

# Path to final .truncspec file
TRUNCSPEC_OUT_PATH = "/home/haig/Repos/micro16s/extract16s/asvs2truncspec_out/new.truncspec"

# Reference sequence IDs (must match Step 1)
ARC_REF_SEQ_ID = "RS_GCF_022846175.1~NZ_AP025587.1-#2"
BAC_REF_SEQ_ID = "RS_GCF_030545895.1~NZ_JAUOMX010000042.1"

# Intermediate directory structure (derived from INFO_OUT_DIR)
INTER_DIR = INFO_OUT_DIR + "/intermediates"
STAGED_DIR = INTER_DIR + "/01_staged_inputs"
HMMER_OUT_DIR = INTER_DIR + "/02_hmmer"
RESULTS_DIR = INTER_DIR + "/03_results"

# Final region plots
FINAL_REGIONS_ARC_PLOT_PATH = INFO_OUT_DIR + "/final_regions_arc.png"
FINAL_REGIONS_BAC_PLOT_PATH = INFO_OUT_DIR + "/final_regions_bac.png"

# Hit filtering thresholds
ARC_HMMER_MAX_EVALUE = 1e-4
BAC_HMMER_MAX_EVALUE = 1e-5
MIN_HIT_COVERAGE_FRAC = 0.6

# Consensus calling parameters
ARC_MIN_HITS_PER_DATASET = 5
BAC_MIN_HITS_PER_DATASET = 5
OUTLIER_TOL_BP = 30

# Region length buffers
MIN_LEN_BUFFER = 50
MAX_LEN_BUFFER = 50

# Optional final processing step 1: Cross‑domain bootstrapping
#    In cases where one domain yeilds coordinates, but the other does not, we can use the 
#    coordinates from the domain that did yeild coordinates to bootstrap the coordinates 
#    for the domain that did not yeild coordinates. It refers to other regions that did 
#    yeild coordinates for both domains. It will only consider relying on other regions 
#    that are within CROSS_DOMAIN_BOOTSTRAPPING_MAX_DISTANCE bp of the region that did 
#    yeild coordinates.
USE_CROSS_DOMAIN_BOOTSTRAPPING = True
CROSS_DOMAIN_BOOTSTRAPPING_MAX_DISTANCE = 90 # bp

# Optional final processing step 2: Region redundancy minimisation
#    If we are using a lot of datasets, we may have a lot of regions that are very similar 
#    and therefore are redundant. This feature merges regions that have both their start 
#    and end coordinates within a certain distance (REGION_REDUNDANCY_MINIMISATION_MAX_DISTANCE) 
#    of each other for *both domains*. 
USE_REGION_REDUNDANCY_MINIMISATION = True
REGION_REDUNDANCY_MINIMISATION_MAX_DISTANCE = 40 # bp
REGION_REDUNDANCY_MINIMISATION_NAME_JOINER = "+"

# Optional final processing step 3: Region renaming
#    Sometimes we might want to rename regions. For example, when using region redundancy minimisation,
#    the names get concatenated resulting in long complicated names. This feature renames the regions.
#    The method "uniform" will rename the regions to a uniform name like "Reg-001", "Reg-002", etc.
#    The method "dataset" will rename the regions based on only the first dataset name in the region name like "{dataset}-001", "{dataset}-002", etc.
#    A helpful TSV file is written which maps the old name to the new name called "region_renamings.tsv".
USE_REGION_RENAMING = True
REGION_RENAMING_METHOD = "uniform" # "uniform" or "dataset"




# =============================================================================
# STOCKHOLM PARSING - Build model→reference mappings
# =============================================================================

def parse_stockholm_for_mapping(sto_path):
    """Parse a Stockholm alignment file to build model_pos → ref_pos mapping.
    
    Walks the alignment to track:
    - model_pos: increments only on match-state columns (RF='x')
    - ref_pos: increments only when ref sequence has a nucleotide (not gap)
    
    Returns dict mapping model_pos (1-based) → ref_pos (1-based).
    Missing mappings (ref has gap at match-state column) are not included.
    """
    # read entire file
    with open(sto_path, 'r') as f:
        lines = f.readlines()
    
    # collect all sequence and RF lines (handle multi-block Stockholm)
    ref_seq_parts = []
    rf_parts = []
    ref_name = None
    
    for line in lines:
        line = line.rstrip('\n')
        
        # skip comments, empty lines, and metadata
        if not line or line.startswith('#=GS') or line.startswith('#=GR') or line.startswith('#=GC PP'):
            continue
        
        # RF annotation line
        if line.startswith('#=GC RF'):
            rf_part = line.split()[-1]  # last token is the RF string
            rf_parts.append(rf_part)
            continue
        
        # skip other #= lines and format markers
        if line.startswith('#') or line == '//' or line.startswith('//'):
            continue
        
        # sequence line: "seqname    ACGU..."
        parts = line.split()
        if len(parts) >= 2:
            seq_name = parts[0]
            seq_part = parts[-1]  # last token is the sequence
            
            # capture the reference sequence name on first encounter
            if ref_name is None:
                ref_name = seq_name
            
            # only process our reference sequence
            if seq_name == ref_name:
                ref_seq_parts.append(seq_part)
    
    # combine all blocks
    ref_seq = ''.join(ref_seq_parts)
    rf_line = ''.join(rf_parts)
    
    if len(ref_seq) != len(rf_line):
        raise ValueError(
            f"Reference sequence length ({len(ref_seq)}) != RF line length ({len(rf_line)}) "
            f"in {sto_path}"
        )
    
    # build model→ref mapping
    mapping = {}
    model_pos = 0  # will be 1-based after increment
    ref_pos = 0    # will be 1-based after increment
    
    for i in range(len(rf_line)):
        rf_char = rf_line[i]
        seq_char = ref_seq[i].upper()
        
        # check if this is a match-state column (RF='x')
        is_match_state = (rf_char == 'x')
        
        # check if ref has a nucleotide (not gap)
        is_nucleotide = seq_char in 'ACGTURYKMSWBDHVN'
        
        # increment ref_pos if nucleotide present
        if is_nucleotide:
            ref_pos += 1
        
        # increment model_pos if match-state column
        if is_match_state:
            model_pos += 1
            # record mapping only if ref has nucleotide at this position
            if is_nucleotide:
                mapping[model_pos] = ref_pos
    
    return mapping


def get_ref_length_from_stockholm(sto_path):
    """Return ungapped reference length from a Stockholm alignment."""
    with open(sto_path, 'r') as f:
        lines = f.readlines()
    
    ref_seq_parts = []
    ref_name = None
    
    for line in lines:
        line = line.rstrip('\n')
        
        if not line or line.startswith('#=GS') or line.startswith('#=GR') or line.startswith('#=GC PP'):
            continue
        
        if line.startswith('#=GC RF'):
            continue
        
        if line.startswith('#') or line == '//' or line.startswith('//'):
            continue
        
        parts = line.split()
        if len(parts) >= 2:
            seq_name = parts[0]
            seq_part = parts[-1]
            
            if ref_name is None:
                ref_name = seq_name
            
            if seq_name == ref_name:
                ref_seq_parts.append(seq_part)
    
    ref_seq = ''.join(ref_seq_parts).upper()
    return sum(1 for c in ref_seq if c in 'ACGTURYKMSWBDHVN')


def snap_to_mapped_position(model_pos, mapping, direction, max_model_pos):
    """Find nearest mapped model position by snapping in a direction.
    
    direction: 'forward' (for start) or 'backward' (for end)
    Returns ref_pos or None if no mapping found.
    """
    if direction == 'forward':
        # search forward from model_pos
        for pos in range(model_pos, max_model_pos + 1):
            if pos in mapping:
                return mapping[pos]
    else:  # backward
        # search backward from model_pos
        for pos in range(model_pos, 0, -1):
            if pos in mapping:
                return mapping[pos]
    
    return None


# =============================================================================
# TBLOUT PARSING - Extract and filter HMMER hits
# =============================================================================

def parse_tblout(tblout_path, max_evalue, min_coverage_frac):
    """Parse a gzipped nhmmer tblout file and filter hits.
    
    Columns in nhmmer tblout (space-separated):
    0: target name
    1: accession
    2: query name
    3: accession
    4: hmmfrom
    5: hmm to
    6: alifrom
    7: ali to
    8: envfrom
    9: env to
    10: sq len
    11: strand
    12: E-value
    13: score
    14: bias
    15+: description
    
    Returns list of dicts with filtered hits:
    {
        'target': full header string,
        'dataset': dataset_id parsed from header,
        'asv': asv_id parsed from header,
        'hmmfrom': int,
        'hmmto': int,
        'envfrom': int,
        'envto': int,
        'sq_len': int,
        'strand': '+' or '-',
        'evalue': float,
        'score': float,
        'coverage': float
    }
    """
    hits = []
    total_hits = 0
    evalue_passed = 0
    coverage_passed = 0
    
    # open gzipped file
    open_func = gzip.open if tblout_path.endswith('.gz') else open
    mode = 'rt' if tblout_path.endswith('.gz') else 'r'
    
    with open_func(tblout_path, mode) as f:
        for line in f:
            # skip comments
            if line.startswith('#'):
                continue
            
            total_hits += 1
            
            # split on whitespace (target name may contain description after)
            parts = line.split()
            if len(parts) < 15:
                continue
            
            # parse fields
            target = parts[0]
            hmmfrom = int(parts[4])
            hmmto = int(parts[5])
            alifrom = int(parts[6])
            alito = int(parts[7])
            envfrom = int(parts[8])
            envto = int(parts[9])
            sq_len = int(parts[10])
            strand = parts[11]
            evalue = float(parts[12])
            score = float(parts[13])
            
            # filter by e-value
            if evalue > max_evalue:
                continue
            evalue_passed += 1
            
            # compute coverage using env coordinates (handle strand)
            env_span = abs(envto - envfrom) + 1
            coverage = env_span / sq_len
            
            # filter by coverage
            if coverage < min_coverage_frac:
                continue
            coverage_passed += 1
            
            # parse dataset and asv from target header
            # format: dataset=XXX|asv=YYY|src=ZZZ.fna
            dataset_id = None
            asv_id = None
            for field in target.split('|'):
                if field.startswith('dataset='):
                    dataset_id = field[8:]
                elif field.startswith('asv='):
                    asv_id = field[4:]
            
            if dataset_id is None or asv_id is None:
                continue  # malformed header
            
            hits.append({
                'target': target,
                'dataset': dataset_id,
                'asv': asv_id,
                'hmmfrom': hmmfrom,
                'hmmto': hmmto,
                'envfrom': envfrom,
                'envto': envto,
                'sq_len': sq_len,
                'strand': strand,
                'evalue': evalue,
                'score': score,
                'coverage': coverage
            })
    
    filter_stats = {
        'total_hits': total_hits,
        'evalue_passed': evalue_passed,
        'coverage_passed': coverage_passed
    }
    
    return hits, filter_stats


def select_best_hit_per_asv(hits):
    """Select the single best hit per ASV (highest bitscore).
    
    Returns dict mapping asv_id → best hit dict.
    Also returns list of secondary hits that were discarded.
    """
    best_hits = {}
    secondary_hits = []
    
    for hit in hits:
        asv = hit['asv']
        
        if asv not in best_hits:
            best_hits[asv] = hit
        else:
            # compare by score (higher is better)
            if hit['score'] > best_hits[asv]['score']:
                secondary_hits.append(best_hits[asv])
                best_hits[asv] = hit
            else:
                secondary_hits.append(hit)
    
    return best_hits, secondary_hits


def apply_best_hmm_gate(arc_best_hits, bac_best_hits):
    """Apply best-HMM gate: keep only the HMM with higher bitscore per ASV.
    
    If an ASV has passing hits to both HMMs, keep only the HMM with higher score.
    Tie-breaker: lower E-value wins; if still tied, prefer bacterial.
    
    Returns tuple: (arc_gated, bac_gated) - dicts mapping asv → hit
    """
    # get all unique ASVs
    all_asvs = set(arc_best_hits.keys()) | set(bac_best_hits.keys())
    
    arc_gated = {}
    bac_gated = {}
    
    for asv in all_asvs:
        arc_hit = arc_best_hits.get(asv)
        bac_hit = bac_best_hits.get(asv)
        
        if arc_hit is None and bac_hit is None:
            continue
        
        if arc_hit is None:
            bac_gated[asv] = bac_hit
            continue
        
        if bac_hit is None:
            arc_gated[asv] = arc_hit
            continue
        
        # both have hits - compare scores
        if arc_hit['score'] > bac_hit['score']:
            arc_gated[asv] = arc_hit
        elif bac_hit['score'] > arc_hit['score']:
            bac_gated[asv] = bac_hit
        else:
            # tie on score - use e-value (lower is better)
            if arc_hit['evalue'] < bac_hit['evalue']:
                arc_gated[asv] = arc_hit
            elif bac_hit['evalue'] < arc_hit['evalue']:
                bac_gated[asv] = bac_hit
            else:
                # still tied - prefer bacterial
                bac_gated[asv] = bac_hit
    
    return arc_gated, bac_gated


# =============================================================================
# COORDINATE CONVERSION - HMM coords to reference coords
# =============================================================================

def convert_hit_to_ref_coords(hit, mapping, max_model_pos):
    """Convert a hit's HMM coordinates to reference coordinates.
    
    Uses snapping to handle gaps in the reference at model positions.
    Returns (ref_start, ref_end) or (None, None) if snapping fails.
    """
    hmmfrom = hit['hmmfrom']
    hmmto = hit['hmmto']
    
    # snap hmmfrom forward to find ref_start
    ref_start = snap_to_mapped_position(hmmfrom, mapping, 'forward', max_model_pos)
    
    # snap hmmto backward to find ref_end
    ref_end = snap_to_mapped_position(hmmto, mapping, 'backward', max_model_pos)
    
    return ref_start, ref_end


def add_ref_coords_to_hits(hits_dict, mapping):
    """Add ref_start and ref_end to each hit in the dict.
    
    Modifies hits in place.
    Returns count of hits with valid ref coords.
    """
    max_model_pos = max(mapping.keys()) if mapping else 0
    valid_count = 0
    
    for asv, hit in hits_dict.items():
        ref_start, ref_end = convert_hit_to_ref_coords(hit, mapping, max_model_pos)
        hit['ref_start'] = ref_start
        hit['ref_end'] = ref_end
        if ref_start is not None and ref_end is not None:
            valid_count += 1
    
    return valid_count


# =============================================================================
# CONSENSUS CALLING - Compute robust consensus regions per dataset
# =============================================================================

def compute_median(values):
    """Compute median. For even length, returns lower of two middle values."""
    if not values:
        return None
    sorted_vals = sorted(values)
    n = len(sorted_vals)
    mid = n // 2
    if n % 2 == 1:
        return sorted_vals[mid]
    else:
        # even length - return lower of two middle values
        return sorted_vals[mid - 1]


def compute_consensus_region(hits, outlier_tol_bp, min_hits):
    """Compute consensus start/end for a list of hits.
    
    Algorithm:
    1. Compute median start and end
    2. Filter out outliers (abs(val - median) > outlier_tol_bp)
    3. Recompute median on remaining hits
    4. Require at least min_hits after filtering
    
    Returns dict with:
    - 'start': consensus start (or None)
    - 'end': consensus end (or None)
    - 'n_hits_initial': number of hits before outlier filtering
    - 'n_hits_final': number of hits after outlier filtering
    - 'hits_used': list of ASV IDs used in final consensus
    """
    # filter to hits with valid ref coords
    valid_hits = [h for h in hits if h['ref_start'] is not None and h['ref_end'] is not None]
    
    if len(valid_hits) < min_hits:
        return {
            'start': None,
            'end': None,
            'n_hits_initial': len(valid_hits),
            'n_hits_final': 0,
            'hits_used': []
        }
    
    # extract start and end values
    starts = [h['ref_start'] for h in valid_hits]
    ends = [h['ref_end'] for h in valid_hits]
    
    # compute initial medians
    start_med = compute_median(starts)
    end_med = compute_median(ends)
    
    # filter outliers
    inlier_hits = []
    for h in valid_hits:
        start_ok = abs(h['ref_start'] - start_med) <= outlier_tol_bp
        end_ok = abs(h['ref_end'] - end_med) <= outlier_tol_bp
        if start_ok and end_ok:
            inlier_hits.append(h)
    
    # check minimum hits after filtering
    if len(inlier_hits) < min_hits:
        return {
            'start': None,
            'end': None,
            'n_hits_initial': len(valid_hits),
            'n_hits_final': len(inlier_hits),
            'hits_used': []
        }
    
    # recompute medians on inliers
    final_starts = [h['ref_start'] for h in inlier_hits]
    final_ends = [h['ref_end'] for h in inlier_hits]
    final_start = compute_median(final_starts)
    final_end = compute_median(final_ends)
    
    return {
        'start': final_start,
        'end': final_end,
        'n_hits_initial': len(valid_hits),
        'n_hits_final': len(inlier_hits),
        'hits_used': [h['asv'] for h in inlier_hits]
    }


def group_hits_by_dataset(hits_dict):
    """Group hits by dataset_id.
    
    Returns dict mapping dataset_id → list of hits.
    """
    by_dataset = {}
    for asv, hit in hits_dict.items():
        ds = hit['dataset']
        if ds not in by_dataset:
            by_dataset[ds] = []
        by_dataset[ds].append(hit)
    return by_dataset


def compute_min_max_len(arc_start, arc_end, bac_start, bac_end, min_len_buffer, max_len_buffer):
    """Compute min_len and max_len based on available domain lengths."""
    lengths = []
    if arc_start is not None and arc_end is not None:
        lengths.append(arc_end - arc_start + 1)
    if bac_start is not None and bac_end is not None:
        lengths.append(bac_end - bac_start + 1)
    
    if not lengths:
        return None, None
    
    min_len = min(lengths) - min_len_buffer
    max_len = max(lengths) + max_len_buffer
    if min_len < 1:
        min_len = 1
    
    return min_len, max_len


def count_complete_regions(dataset_calls):
    """Count datasets with both arc and bac coords present."""
    return sum(
        1 for call in dataset_calls.values()
        if call['arc_start'] is not None and call['arc_end'] is not None
        and call['bac_start'] is not None and call['bac_end'] is not None
    )


def build_bootstrap_reference_regions(dataset_calls):
    """Collect datasets that have complete regions for both domains."""
    refs = []
    for ds_id, call in dataset_calls.items():
        if call['arc_start'] is None or call['arc_end'] is None:
            continue
        if call['bac_start'] is None or call['bac_end'] is None:
            continue
        refs.append({
            'dataset_id': ds_id,
            'arc_start': call['arc_start'],
            'arc_end': call['arc_end'],
            'bac_start': call['bac_start'],
            'bac_end': call['bac_end']
        })
    return refs


def find_closest_reference(ref_regions, domain, coord_type, target_coord):
    """Find closest reference region by coord for a given domain and coord type."""
    best = None
    best_abs = None
    
    key = f"{domain}_{coord_type}"
    for ref in ref_regions:
        ref_coord = ref.get(key)
        if ref_coord is None:
            continue
        dist = target_coord - ref_coord
        abs_dist = abs(dist)
        if best_abs is None or abs_dist < best_abs:
            best = (ref, dist, abs_dist)
            best_abs = abs_dist
        elif abs_dist == best_abs:
            if ref['dataset_id'] < best[0]['dataset_id']:
                best = (ref, dist, abs_dist)
    
    return best


def apply_cross_domain_bootstrapping(dataset_calls, max_distance):
    """Fill missing domain coords by referencing nearby complete regions."""
    ref_regions = build_bootstrap_reference_regions(dataset_calls)
    
    log_entries = []
    bootstrapped_starts = 0
    bootstrapped_ends = 0
    
    for ds_id, call in dataset_calls.items():
        arc_missing = call['arc_start'] is None or call['arc_end'] is None
        bac_missing = call['bac_start'] is None or call['bac_end'] is None
        
        if arc_missing and bac_missing:
            continue
        
        if arc_missing and not bac_missing:
            missing_domain = 'arc'
            known_domain = 'bac'
        elif bac_missing and not arc_missing:
            missing_domain = 'bac'
            known_domain = 'arc'
        else:
            continue
        
        for coord_type in ['start', 'end']:
            missing_key = f"{missing_domain}_{coord_type}"
            known_key = f"{known_domain}_{coord_type}"
            
            if call[missing_key] is not None:
                continue
            
            known_coord = call.get(known_key)
            if known_coord is None:
                log_entries.append({
                    'dataset_id': ds_id,
                    'missing_domain': missing_domain,
                    'coord_type': coord_type,
                    'known_coord': 'NA',
                    'ref_dataset': 'NA',
                    'ref_known_coord': 'NA',
                    'ref_missing_coord': 'NA',
                    'distance': 'NA',
                    'bootstrapped_coord': 'NA',
                    'status': 'failed',
                    'note': 'known_coord_missing'
                })
                continue
            
            closest = find_closest_reference(ref_regions, known_domain, coord_type, known_coord)
            if closest is None:
                log_entries.append({
                    'dataset_id': ds_id,
                    'missing_domain': missing_domain,
                    'coord_type': coord_type,
                    'known_coord': known_coord,
                    'ref_dataset': 'NA',
                    'ref_known_coord': 'NA',
                    'ref_missing_coord': 'NA',
                    'distance': 'NA',
                    'bootstrapped_coord': 'NA',
                    'status': 'failed',
                    'note': 'no_reference_regions'
                })
                continue
            
            ref_region, dist, abs_dist = closest
            if abs_dist > max_distance:
                log_entries.append({
                    'dataset_id': ds_id,
                    'missing_domain': missing_domain,
                    'coord_type': coord_type,
                    'known_coord': known_coord,
                    'ref_dataset': ref_region['dataset_id'],
                    'ref_known_coord': ref_region[f"{known_domain}_{coord_type}"],
                    'ref_missing_coord': ref_region[f"{missing_domain}_{coord_type}"],
                    'distance': dist,
                    'bootstrapped_coord': 'NA',
                    'status': 'failed',
                    'note': 'distance_exceeded'
                })
                continue
            
            ref_missing_coord = ref_region.get(f"{missing_domain}_{coord_type}")
            if ref_missing_coord is None:
                log_entries.append({
                    'dataset_id': ds_id,
                    'missing_domain': missing_domain,
                    'coord_type': coord_type,
                    'known_coord': known_coord,
                    'ref_dataset': ref_region['dataset_id'],
                    'ref_known_coord': ref_region[f"{known_domain}_{coord_type}"],
                    'ref_missing_coord': 'NA',
                    'distance': dist,
                    'bootstrapped_coord': 'NA',
                    'status': 'failed',
                    'note': 'reference_missing_coord'
                })
                continue
            
            bootstrapped_coord = ref_missing_coord + dist
            if bootstrapped_coord < 1:
                log_entries.append({
                    'dataset_id': ds_id,
                    'missing_domain': missing_domain,
                    'coord_type': coord_type,
                    'known_coord': known_coord,
                    'ref_dataset': ref_region['dataset_id'],
                    'ref_known_coord': ref_region[f"{known_domain}_{coord_type}"],
                    'ref_missing_coord': ref_missing_coord,
                    'distance': dist,
                    'bootstrapped_coord': bootstrapped_coord,
                    'status': 'failed',
                    'note': 'bootstrapped_coord_out_of_range'
                })
                continue
            
            call[missing_key] = bootstrapped_coord
            if 'bootstrapped_coords' not in call:
                call['bootstrapped_coords'] = set()
            call['bootstrapped_coords'].add(missing_key)
            
            log_entries.append({
                'dataset_id': ds_id,
                'missing_domain': missing_domain,
                'coord_type': coord_type,
                'known_coord': known_coord,
                'ref_dataset': ref_region['dataset_id'],
                'ref_known_coord': ref_region[f"{known_domain}_{coord_type}"],
                'ref_missing_coord': ref_missing_coord,
                'distance': dist,
                'bootstrapped_coord': bootstrapped_coord,
                'status': 'bootstrapped',
                'note': ''
            })
            
            if coord_type == 'start':
                bootstrapped_starts += 1
            else:
                bootstrapped_ends += 1

    # sanity check: ensure start <= end for any bootstrapped domain
    for ds_id, call in dataset_calls.items():
        for domain in ['arc', 'bac']:
            start_key = f"{domain}_start"
            end_key = f"{domain}_end"
            start = call.get(start_key)
            end = call.get(end_key)
            if start is not None and end is not None and start > end:
                call[start_key] = None
                call[end_key] = None
                if 'bootstrapped_coords' in call:
                    call['bootstrapped_coords'].discard(start_key)
                    call['bootstrapped_coords'].discard(end_key)
                log_entries.append({
                    'dataset_id': ds_id,
                    'missing_domain': domain,
                    'coord_type': 'both',
                    'known_coord': 'NA',
                    'ref_dataset': 'NA',
                    'ref_known_coord': 'NA',
                    'ref_missing_coord': 'NA',
                    'distance': 'NA',
                    'bootstrapped_coord': 'NA',
                    'status': 'failed',
                    'note': 'start_gt_end_after_bootstrapping'
                })
    
    return dataset_calls, log_entries, bootstrapped_starts, bootstrapped_ends


# =============================================================================
# REGION REDUNDANCY MINIMISATION
# =============================================================================

def round_mean(values):
    """Return the mean of values rounded to nearest int (0.5 rounds up)."""
    return int((sum(values) / len(values)) + 0.5)


def regions_within_distance(call_a, call_b, max_distance):
    """Check if two complete regions are within max_distance for both domains."""
    for key in ['arc_start', 'arc_end', 'bac_start', 'bac_end']:
        if abs(call_a[key] - call_b[key]) > max_distance:
            return False
    return True


def region_within_group(candidate, group, max_distance):
    """Require candidate to be within distance of every member in the group."""
    for member in group:
        if not regions_within_distance(candidate, member, max_distance):
            return False
    return True


def apply_region_redundancy_minimisation(dataset_calls, max_distance, name_joiner):
    """Merge similar regions based on arc/bac start/end proximity."""
    complete_entries = []
    incomplete_calls = {}
    
    for ds_id, call in dataset_calls.items():
        if (call['arc_start'] is None or call['arc_end'] is None or
                call['bac_start'] is None or call['bac_end'] is None):
            incomplete_calls[ds_id] = copy.deepcopy(call)
            continue
        
        entry = {
            'dataset_id': ds_id,
            'arc_start': call['arc_start'],
            'arc_end': call['arc_end'],
            'bac_start': call['bac_start'],
            'bac_end': call['bac_end']
        }
        complete_entries.append(entry)
    
    remaining = sorted(
        complete_entries,
        key=lambda e: (e['arc_start'], e['arc_end'], e['bac_start'], e['bac_end'], e['dataset_id'])
    )
    
    merged_calls = {}
    merge_groups = []
    
    while remaining:
        seed = remaining.pop(0)
        group = [seed]
        
        changed = True
        while changed:
            changed = False
            for entry in list(remaining):
                if region_within_group(entry, group, max_distance):
                    group.append(entry)
                    remaining.remove(entry)
                    changed = True
        
        if len(group) == 1:
            ds_id = seed['dataset_id']
            merged_calls[ds_id] = copy.deepcopy(dataset_calls[ds_id])
            continue
        
        group_ids = sorted(entry['dataset_id'] for entry in group)
        merged_id = name_joiner.join(group_ids)
        
        arc_start = round_mean([entry['arc_start'] for entry in group])
        arc_end = round_mean([entry['arc_end'] for entry in group])
        bac_start = round_mean([entry['bac_start'] for entry in group])
        bac_end = round_mean([entry['bac_end'] for entry in group])
        
        min_len, max_len = compute_min_max_len(
            arc_start, arc_end, bac_start, bac_end, MIN_LEN_BUFFER, MAX_LEN_BUFFER
        )
        
        merged_calls[merged_id] = {
            'arc_start': arc_start,
            'arc_end': arc_end,
            'bac_start': bac_start,
            'bac_end': bac_end,
            'min_len': min_len,
            'max_len': max_len
        }
        
        merge_groups.append(group_ids)
    
    for ds_id, call in incomplete_calls.items():
        merged_calls[ds_id] = copy.deepcopy(call)
    
    return merged_calls, merge_groups


# =============================================================================
# REGION RENAMING
# =============================================================================

def select_dataset_prefix(region_name, name_joiner):
    """Select dataset prefix from region name for dataset-based renaming."""
    parts = region_name.split(name_joiner)
    
    non_hmc_parts = [p for p in parts if 'HMC' not in p]
    
    if non_hmc_parts:
        return non_hmc_parts[0]
    else:
        return parts[0]


def apply_region_renaming(dataset_calls, method, name_joiner):
    """Rename regions according to the specified method.
    
    Returns:
    - renamed_calls: dict with new names as keys
    - renaming_map: dict mapping old_name -> new_name
    """
    sorted_names = sorted(dataset_calls.keys())
    
    renaming_map = {}
    renamed_calls = {}
    
    global_counter = 1
    
    for old_name in sorted_names:
        if method == 'uniform':
            new_name = f"Reg-{global_counter:03d}"
        elif method == 'dataset':
            prefix = select_dataset_prefix(old_name, name_joiner)
            new_name = f"{prefix}-{global_counter:03d}"
        else:
            raise ValueError(f"Unknown renaming method: {method}")
        
        renaming_map[old_name] = new_name
        renamed_calls[new_name] = copy.deepcopy(dataset_calls[old_name])
        global_counter += 1
    
    return renamed_calls, renaming_map


def write_renaming_map(renaming_map, out_path):
    """Write region renaming map to TSV file."""
    with open(out_path, 'w') as f:
        f.write("old_name\tnew_name\n")
        for old_name in sorted(renaming_map.keys()):
            f.write(f"{old_name}\t{renaming_map[old_name]}\n")


# =============================================================================
# OUTPUT WRITING
# =============================================================================

def insert_suffix_before_ext(path, suffix):
    """Insert suffix before extension (or append if no extension)."""
    base, ext = os.path.splitext(path)
    if not ext:
        return base + suffix
    return base + suffix + ext


def write_bootstrap_log(log_entries, out_path):
    """Write bootstrapping decisions to TSV file."""
    with open(out_path, 'w') as f:
        f.write("dataset_id\tmissing_domain\tcoord_type\tknown_coord\tref_dataset\t")
        f.write("ref_known_coord\tref_missing_coord\tdistance\tbootstrapped_coord\tstatus\tnote\n")
        for entry in log_entries:
            f.write(
                f"{entry['dataset_id']}\t{entry['missing_domain']}\t{entry['coord_type']}\t"
                f"{entry['known_coord']}\t{entry['ref_dataset']}\t{entry['ref_known_coord']}\t"
                f"{entry['ref_missing_coord']}\t{entry['distance']}\t{entry['bootstrapped_coord']}\t"
                f"{entry['status']}\t{entry['note']}\n"
            )

def write_model_to_ref_mapping(mapping, out_path):
    """Write model→ref mapping to TSV file."""
    with open(out_path, 'w') as f:
        f.write("model_pos\tref_pos\n")
        for model_pos in sorted(mapping.keys()):
            f.write(f"{model_pos}\t{mapping[model_pos]}\n")


def write_asv_hits_tsv(hits_dict, out_path, domain):
    """Write per-ASV hits to TSV file."""
    with open(out_path, 'w') as f:
        f.write("asv_id\tdataset_id\tdomain\thmmfrom\thmmto\tref_start\tref_end\t")
        f.write("envfrom\tenvto\tsq_len\tstrand\tevalue\tscore\tcoverage\n")
        
        for asv in sorted(hits_dict.keys()):
            hit = hits_dict[asv]
            ref_start = hit.get('ref_start', 'NA')
            ref_end = hit.get('ref_end', 'NA')
            if ref_start is None:
                ref_start = 'NA'
            if ref_end is None:
                ref_end = 'NA'
            
            f.write(f"{asv}\t{hit['dataset']}\t{domain}\t")
            f.write(f"{hit['hmmfrom']}\t{hit['hmmto']}\t{ref_start}\t{ref_end}\t")
            f.write(f"{hit['envfrom']}\t{hit['envto']}\t{hit['sq_len']}\t{hit['strand']}\t")
            f.write(f"{hit['evalue']:.2e}\t{hit['score']:.1f}\t{hit['coverage']:.3f}\n")


def write_dataset_region_calls(dataset_calls, out_path):
    """Write dataset region calls to TSV file."""
    with open(out_path, 'w') as f:
        f.write("dataset_id\tarc_start\tarc_end\tbac_start\tbac_end\tmin_len\tmax_len\n")
        
        for ds_id in sorted(dataset_calls.keys()):
            call = dataset_calls[ds_id]
            arc_start = call['arc_start'] if call['arc_start'] is not None else 'NA'
            arc_end = call['arc_end'] if call['arc_end'] is not None else 'NA'
            bac_start = call['bac_start'] if call['bac_start'] is not None else 'NA'
            bac_end = call['bac_end'] if call['bac_end'] is not None else 'NA'
            min_len = call['min_len'] if call['min_len'] is not None else 'NA'
            max_len = call['max_len'] if call['max_len'] is not None else 'NA'
            
            f.write(f"{ds_id}\t{arc_start}\t{arc_end}\t{bac_start}\t{bac_end}\t{min_len}\t{max_len}\n")


def write_dataset_region_qc(dataset_qc, out_path):
    """Write detailed QC info per dataset to TSV file."""
    with open(out_path, 'w') as f:
        f.write("dataset_id\tdomain\tn_hits_initial\tn_hits_final\tstart\tend\n")
        
        for ds_id in sorted(dataset_qc.keys()):
            for domain in ['arc', 'bac']:
                info = dataset_qc[ds_id].get(domain, {})
                n_init = info.get('n_hits_initial', 0)
                n_final = info.get('n_hits_final', 0)
                start = info.get('start', 'NA')
                end = info.get('end', 'NA')
                if start is None:
                    start = 'NA'
                if end is None:
                    end = 'NA'
                
                f.write(f"{ds_id}\t{domain}\t{n_init}\t{n_final}\t{start}\t{end}\n")


def write_filter_stats(arc_stats, bac_stats, arc_gated, bac_gated, total_asvs, out_path):
    """Write filter statistics to TSV file."""
    # count strand distribution
    arc_plus = sum(1 for h in arc_gated.values() if h['strand'] == '+')
    arc_minus = len(arc_gated) - arc_plus
    bac_plus = sum(1 for h in bac_gated.values() if h['strand'] == '+')
    bac_minus = len(bac_gated) - bac_plus
    
    with open(out_path, 'w') as f:
        f.write("statistic\tarc\tbac\n")
        f.write(f"total_asvs_searched\t{total_asvs}\t{total_asvs}\n")
        f.write(f"total_hits_reported\t{arc_stats['total_hits']}\t{bac_stats['total_hits']}\n")
        f.write(f"hits_passing_evalue\t{arc_stats['evalue_passed']}\t{bac_stats['evalue_passed']}\n")
        f.write(f"hits_passing_evalue_and_coverage\t{arc_stats['coverage_passed']}\t{bac_stats['coverage_passed']}\n")
        f.write(f"asvs_with_hit_after_best_hit_selection\t{len(arc_gated) + len(bac_gated) - len(set(arc_gated.keys()) & set(bac_gated.keys()))}\t-\n")
        f.write(f"asvs_retained_after_best_hmm_gate\t{len(arc_gated)}\t{len(bac_gated)}\n")
        f.write(f"strand_plus\t{arc_plus}\t{bac_plus}\n")
        f.write(f"strand_minus\t{arc_minus}\t{bac_minus}\n")
        
        # compute percent of ASVs with passing hit
        arc_pct = 100 * len(arc_gated) / total_asvs if total_asvs > 0 else 0
        bac_pct = 100 * len(bac_gated) / total_asvs if total_asvs > 0 else 0
        f.write(f"pct_asvs_with_passing_hit\t{arc_pct:.2f}\t{bac_pct:.2f}\n")


def truncate_region_name(region_name):
    """Shorten long region names for plotting."""
    if len(region_name) > 25:
        return region_name[:22] + "..."
    return region_name


def get_complete_regions(dataset_calls):
    """Collect regions that have complete coordinates for both domains."""
    regions = []
    for ds_id in sorted(dataset_calls.keys()):
        call = dataset_calls[ds_id]
        arc_ok = call['arc_start'] is not None and call['arc_end'] is not None
        bac_ok = call['bac_start'] is not None and call['bac_end'] is not None
        if arc_ok and bac_ok:
            regions.append({
                'dataset_id': ds_id,
                'arc_start': call['arc_start'],
                'arc_end': call['arc_end'],
                'bac_start': call['bac_start'],
                'bac_end': call['bac_end']
            })
    return regions


def plot_final_regions(complete_regions, domain, ref_length, out_path):
    """Plot final complete regions across a reference sequence."""
    if not complete_regions:
        return
    
    fig_height = max(2.5, 0.25 * len(complete_regions) + 1.0)
    fig, ax = plt.subplots(figsize=(20, fig_height))
    
    y_positions = list(range(len(complete_regions)))
    labels = [truncate_region_name(r['dataset_id']) for r in complete_regions]
    
    for idx, region in enumerate(complete_regions):
        start = region[f"{domain}_start"]
        end = region[f"{domain}_end"]
        width = end - start + 1
        ax.barh(idx, width, left=start, height=0.8, color="#1f77b4", zorder=2)
    
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Reference position (bp)", fontsize=8)
    ax.set_ylabel("Region", fontsize=8)
    ax.set_xlim(1, ref_length)
    ax.set_xticks(list(range(0, ref_length + 1, 25)))
    ax.tick_params(axis='x', labelsize=7)
    ax.set_axisbelow(True)
    ax.grid(axis='both', linestyle='-', linewidth=0.6, color='#666666', alpha=0.7, zorder=0)
    
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def format_bootstrap_comment(call):
    """Format bootstrapping comment for truncspec lines."""
    coords = call.get('bootstrapped_coords', set())
    if not coords:
        return ''
    
    parts = []
    for domain_key, label in [('arc', 'Archaea'), ('bac', 'Bacteria')]:
        start_key = f"{domain_key}_start"
        end_key = f"{domain_key}_end"
        has_start = start_key in coords
        has_end = end_key in coords
        if not has_start and not has_end:
            continue
        
        if has_start and has_end:
            coord_phrase = "start and end"
        elif has_start:
            coord_phrase = "start"
        else:
            coord_phrase = "end"
        
        parts.append(f"{label} {coord_phrase}")
    
    if not parts:
        return ''
    
    return " # Bootstrapped: " + "; ".join(parts)


def write_truncspec(dataset_calls, arc_ref_id, bac_ref_id, out_path, include_bootstrap_comments=True):
    """Write the final .truncspec file."""
    with open(out_path, 'w') as f:
        # header comments
        f.write("# 16S Region Truncation Specification File (.truncspec)\n")
        f.write("# Generated by asvs2truncspec\n")
        f.write(f"# Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("# This file specifies regions to be extracted from 16S alignments using `extract16s.sh`\n")
        f.write("\n")
        
        # reference sequence IDs
        f.write(f"ARC_REF_SEQ_ID: {arc_ref_id}\n")
        f.write(f"BAC_REF_SEQ_ID: {bac_ref_id}\n")
        
        # region definitions
        for ds_id in sorted(dataset_calls.keys()):
            call = dataset_calls[ds_id]
            
            # check if all coordinates are valid
            arc_ok = call['arc_start'] is not None and call['arc_end'] is not None
            bac_ok = call['bac_start'] is not None and call['bac_end'] is not None
            all_ok = arc_ok and bac_ok
            
            # format values
            arc_start = call['arc_start'] if call['arc_start'] is not None else 'NA'
            arc_end = call['arc_end'] if call['arc_end'] is not None else 'NA'
            bac_start = call['bac_start'] if call['bac_start'] is not None else 'NA'
            bac_end = call['bac_end'] if call['bac_end'] is not None else 'NA'
            min_len = call['min_len'] if call['min_len'] is not None else 'NA'
            max_len = call['max_len'] if call['max_len'] is not None else 'NA'
            
            # format the line
            line = f"{ds_id}: arc_start={arc_start}, arc_end={arc_end}, "
            line += f"bac_start={bac_start}, bac_end={bac_end}, "
            line += f"min_len={min_len}, max_len={max_len}"
            
            # comment out if incomplete
            if not all_ok:
                line = "#" + line
            
            if include_bootstrap_comments:
                line += format_bootstrap_comment(call)
            
            f.write(line + "\n")


def write_about_file(out_path, config_info, stats_info):
    """Write summary about file for Step 3."""
    with open(out_path, 'w') as f:
        f.write("About Processing (Step 3)\n")
        f.write(f"Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write("\n\n")
        
        f.write("Configuration\n")
        f.write("-" * 40 + "\n")
        for key, val in config_info.items():
            f.write(f"{key} = {val}\n")
        f.write("\n\n")
        
        f.write("Processing Summary\n")
        f.write("-" * 40 + "\n")
        for key, val in stats_info.items():
            f.write(f"{key}: {val}\n")


# =============================================================================
# MAIN SCRIPT
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("asvs2truncspec Step 3: Process HMMER outputs")
    print("=" * 70)
    print()
    
    # --- 1. Create output directory ---
    print("Creating output directory...")
    os.makedirs(RESULTS_DIR, exist_ok=True)
    print(f"  Created: {RESULTS_DIR}")
    print()
    
    # --- 2. Parse Stockholm files to build model→ref mappings ---
    print("Parsing Stockholm alignments for model→reference mappings...")
    
    arc_sto_path = HMMER_OUT_DIR + "/arc_ref.sto"
    bac_sto_path = HMMER_OUT_DIR + "/bac_ref.sto"
    
    arc_mapping = parse_stockholm_for_mapping(arc_sto_path)
    print(f"  Archaeal mapping: {len(arc_mapping)} model positions → ref positions")
    
    bac_mapping = parse_stockholm_for_mapping(bac_sto_path)
    print(f"  Bacterial mapping: {len(bac_mapping)} model positions → ref positions")
    
    arc_ref_len = get_ref_length_from_stockholm(arc_sto_path)
    bac_ref_len = get_ref_length_from_stockholm(bac_sto_path)
    
    # write mapping files
    arc_map_path = RESULTS_DIR + "/arc_model_to_ref.tsv"
    bac_map_path = RESULTS_DIR + "/bac_model_to_ref.tsv"
    write_model_to_ref_mapping(arc_mapping, arc_map_path)
    write_model_to_ref_mapping(bac_mapping, bac_map_path)
    print(f"  Wrote: {arc_map_path}")
    print(f"  Wrote: {bac_map_path}")
    print()
    
    # --- 3. Load dataset manifest to get total ASV count and dataset list ---
    print("Loading dataset manifest...")
    manifest_path = STAGED_DIR + "/dataset_manifest.tsv"
    
    dataset_to_asvs = {}  # dataset_id → set of asv_ids
    total_asvs = 0
    
    with open(manifest_path, 'r') as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                ds_id = parts[0]
                asv_id = parts[1]
                if ds_id not in dataset_to_asvs:
                    dataset_to_asvs[ds_id] = set()
                dataset_to_asvs[ds_id].add(asv_id)
                total_asvs += 1
    
    all_datasets = sorted(dataset_to_asvs.keys())
    print(f"  Total ASVs: {total_asvs}")
    print(f"  Total datasets: {len(all_datasets)}")
    print()
    
    # --- 4. Parse tblout files and filter hits ---
    print("Parsing and filtering HMMER tblout files...")
    print(f"  E-value threshold (arc): {ARC_HMMER_MAX_EVALUE}")
    print(f"  E-value threshold (bac): {BAC_HMMER_MAX_EVALUE}")
    print(f"  Min coverage fraction: {MIN_HIT_COVERAGE_FRAC}")
    
    arc_tblout_path = HMMER_OUT_DIR + "/arc.tblout.gz"
    bac_tblout_path = HMMER_OUT_DIR + "/bac.tblout.gz"
    
    arc_hits, arc_filter_stats = parse_tblout(arc_tblout_path, ARC_HMMER_MAX_EVALUE, MIN_HIT_COVERAGE_FRAC)
    print(f"  Archaeal: {arc_filter_stats['total_hits']} total → "
          f"{arc_filter_stats['evalue_passed']} E-value → "
          f"{arc_filter_stats['coverage_passed']} coverage")
    
    bac_hits, bac_filter_stats = parse_tblout(bac_tblout_path, BAC_HMMER_MAX_EVALUE, MIN_HIT_COVERAGE_FRAC)
    print(f"  Bacterial: {bac_filter_stats['total_hits']} total → "
          f"{bac_filter_stats['evalue_passed']} E-value → "
          f"{bac_filter_stats['coverage_passed']} coverage")
    print()
    
    # --- 5. Select best hit per ASV per HMM ---
    print("Selecting best hit per ASV per HMM...")
    
    arc_best, arc_secondary = select_best_hit_per_asv(arc_hits)
    bac_best, bac_secondary = select_best_hit_per_asv(bac_hits)
    
    print(f"  Archaeal: {len(arc_best)} ASVs with best hits ({len(arc_secondary)} secondary discarded)")
    print(f"  Bacterial: {len(bac_best)} ASVs with best hits ({len(bac_secondary)} secondary discarded)")
    print()
    
    # --- 6. Apply best-HMM gate ---
    print("Applying best-HMM gate (keep only winning HMM per ASV)...")
    
    arc_gated, bac_gated = apply_best_hmm_gate(arc_best, bac_best)
    
    print(f"  Archaeal: {len(arc_gated)} ASVs retained")
    print(f"  Bacterial: {len(bac_gated)} ASVs retained")
    print()
    
    # --- 7. Convert HMM coordinates to reference coordinates ---
    print("Converting HMM coordinates to reference coordinates...")
    
    arc_valid = add_ref_coords_to_hits(arc_gated, arc_mapping)
    bac_valid = add_ref_coords_to_hits(bac_gated, bac_mapping)
    
    print(f"  Archaeal: {arc_valid}/{len(arc_gated)} with valid ref coords")
    print(f"  Bacterial: {bac_valid}/{len(bac_gated)} with valid ref coords")
    print()
    
    # --- 8. Write per-ASV hit files ---
    print("Writing per-ASV hit files...")
    
    arc_hits_path = RESULTS_DIR + "/asv_hits_arc.tsv"
    bac_hits_path = RESULTS_DIR + "/asv_hits_bac.tsv"
    
    write_asv_hits_tsv(arc_gated, arc_hits_path, 'arc')
    write_asv_hits_tsv(bac_gated, bac_hits_path, 'bac')
    
    print(f"  Wrote: {arc_hits_path}")
    print(f"  Wrote: {bac_hits_path}")
    print()
    
    # --- 9. Compute consensus regions per dataset ---
    print("Computing consensus regions per dataset...")
    print(f"  Outlier tolerance: {OUTLIER_TOL_BP} bp")
    print(f"  Min hits per dataset (arc): {ARC_MIN_HITS_PER_DATASET}")
    print(f"  Min hits per dataset (bac): {BAC_MIN_HITS_PER_DATASET}")
    
    # group hits by dataset
    arc_by_dataset = group_hits_by_dataset(arc_gated)
    bac_by_dataset = group_hits_by_dataset(bac_gated)
    
    # compute consensus for each dataset
    dataset_calls = {}
    dataset_qc = {}
    
    for ds_id in all_datasets:
        # archaeal consensus
        arc_ds_hits = arc_by_dataset.get(ds_id, [])
        arc_consensus = compute_consensus_region(arc_ds_hits, OUTLIER_TOL_BP, ARC_MIN_HITS_PER_DATASET)
        
        # bacterial consensus
        bac_ds_hits = bac_by_dataset.get(ds_id, [])
        bac_consensus = compute_consensus_region(bac_ds_hits, OUTLIER_TOL_BP, BAC_MIN_HITS_PER_DATASET)
        
        # store QC info
        dataset_qc[ds_id] = {
            'arc': arc_consensus,
            'bac': bac_consensus
        }
        
        # compute min_len and max_len
        arc_start = arc_consensus['start']
        arc_end = arc_consensus['end']
        bac_start = bac_consensus['start']
        bac_end = bac_consensus['end']
        
        min_len, max_len = compute_min_max_len(
            arc_start, arc_end, bac_start, bac_end, MIN_LEN_BUFFER, MAX_LEN_BUFFER
        )
        
        dataset_calls[ds_id] = {
            'arc_start': arc_start,
            'arc_end': arc_end,
            'bac_start': bac_start,
            'bac_end': bac_end,
            'min_len': min_len,
            'max_len': max_len
        }
    
    # count complete vs incomplete before bootstrapping
    complete_count = count_complete_regions(dataset_calls)
    incomplete_count = len(dataset_calls) - complete_count
    
    print(f"  Complete regions (pre-bootstrapping): {complete_count}")
    print(f"  Incomplete regions (pre-bootstrapping): {incomplete_count}")
    print()

    # --- 9b. Optional cross-domain bootstrapping ---
    use_bootstrapping = USE_CROSS_DOMAIN_BOOTSTRAPPING
    boot_log_entries = []
    boot_starts = 0
    boot_ends = 0
    
    if use_bootstrapping and len(all_datasets) < 2:
        print("WARNING: Only one dataset found; disabling cross-domain bootstrapping.")
        use_bootstrapping = False
    
    if use_bootstrapping:
        print("Applying cross-domain bootstrapping...")
        print(f"  Max distance: {CROSS_DOMAIN_BOOTSTRAPPING_MAX_DISTANCE} bp")
        
        dataset_calls_unbootstrapped = copy.deepcopy(dataset_calls)
        
        dataset_calls, boot_log_entries, boot_starts, boot_ends = apply_cross_domain_bootstrapping(
            dataset_calls, CROSS_DOMAIN_BOOTSTRAPPING_MAX_DISTANCE
        )
        
        # recompute min_len/max_len after bootstrapping
        for call in dataset_calls.values():
            call['min_len'], call['max_len'] = compute_min_max_len(
                call['arc_start'], call['arc_end'],
                call['bac_start'], call['bac_end'],
                MIN_LEN_BUFFER, MAX_LEN_BUFFER
            )
        
        complete_count = count_complete_regions(dataset_calls)
        incomplete_count = len(dataset_calls) - complete_count
        
        print(f"  Bootstrapped starts: {boot_starts}")
        print(f"  Bootstrapped ends: {boot_ends}")
        print(f"  Complete regions (post-bootstrapping): {complete_count}")
        print(f"  Incomplete regions (post-bootstrapping): {incomplete_count}")
        print()
    else:
        dataset_calls_unbootstrapped = None
    
    dataset_calls_pre_redundancy = dataset_calls
    complete_count_pre_redundancy = complete_count
    incomplete_count_pre_redundancy = incomplete_count
    
    # --- 10. Write QC and stats files ---
    print("Writing QC and statistics files...")
    
    # dataset region calls
    calls_path = RESULTS_DIR + "/dataset_region_calls.tsv"
    write_dataset_region_calls(dataset_calls_pre_redundancy, calls_path)
    print(f"  Wrote: {calls_path}")
    
    # detailed QC
    qc_path = RESULTS_DIR + "/dataset_region_qc.tsv"
    write_dataset_region_qc(dataset_qc, qc_path)
    print(f"  Wrote: {qc_path}")
    
    # filter stats
    stats_path = RESULTS_DIR + "/hmmer_filter_stats.tsv"
    write_filter_stats(arc_filter_stats, bac_filter_stats, arc_gated, bac_gated, total_asvs, stats_path)
    print(f"  Wrote: {stats_path}")

    # bootstrapping log
    if use_bootstrapping:
        boot_log_path = RESULTS_DIR + "/cross_domain_bootstrapping.tsv"
        write_bootstrap_log(boot_log_entries, boot_log_path)
        print(f"  Wrote: {boot_log_path}")
    print()
    
    # --- 10b. Optional region redundancy minimisation ---
    use_redundancy_minimisation = USE_REGION_REDUNDANCY_MINIMISATION
    redundancy_merge_groups = []
    
    if use_redundancy_minimisation:
        print("Applying region redundancy minimisation...")
        print(f"  Max distance: {REGION_REDUNDANCY_MINIMISATION_MAX_DISTANCE} bp")
        
        dataset_calls_redundancy_min, redundancy_merge_groups = apply_region_redundancy_minimisation(
            dataset_calls_pre_redundancy,
            REGION_REDUNDANCY_MINIMISATION_MAX_DISTANCE,
            REGION_REDUNDANCY_MINIMISATION_NAME_JOINER
        )
        
        complete_count_post_redundancy = count_complete_regions(dataset_calls_redundancy_min)
        incomplete_count_post_redundancy = len(dataset_calls_redundancy_min) - complete_count_post_redundancy
        
        merged_group_count = len(redundancy_merge_groups)
        merged_region_count = sum(len(group) for group in redundancy_merge_groups)
        removed_region_count = merged_region_count - merged_group_count
        
        print(f"  Merge groups: {merged_group_count}")
        print(f"  Regions merged: {merged_region_count} -> {merged_group_count} (removed {removed_region_count})")
        print(f"  Complete regions (post-redundancy): {complete_count_post_redundancy}")
        print(f"  Incomplete regions (post-redundancy): {incomplete_count_post_redundancy}")
        print()
    else:
        dataset_calls_redundancy_min = dataset_calls_pre_redundancy
        complete_count_post_redundancy = complete_count_pre_redundancy
        incomplete_count_post_redundancy = incomplete_count_pre_redundancy
    
    # --- 10c. Optional region renaming ---
    use_region_renaming = USE_REGION_RENAMING
    renaming_map = {}
    
    if use_region_renaming:
        print("Applying region renaming...")
        print(f"  Method: {REGION_RENAMING_METHOD}")
        
        dataset_calls_pre_renaming = dataset_calls_redundancy_min
        
        dataset_calls_renamed, renaming_map = apply_region_renaming(
            dataset_calls_redundancy_min,
            REGION_RENAMING_METHOD,
            REGION_REDUNDANCY_MINIMISATION_NAME_JOINER
        )
        
        print(f"  Renamed {len(renaming_map)} regions")
        print()
    else:
        dataset_calls_renamed = dataset_calls_redundancy_min
        dataset_calls_pre_renaming = None
    
    # --- 11. Write final .truncspec file ---
    print("Writing final .truncspec file...")
    if use_bootstrapping:
        unbootstrapped_path = insert_suffix_before_ext(TRUNCSPEC_OUT_PATH, "_no-bootstrapping")
        write_truncspec(dataset_calls_unbootstrapped, ARC_REF_SEQ_ID, BAC_REF_SEQ_ID, unbootstrapped_path)
        print(f"  Wrote: {unbootstrapped_path}")
    
    if use_redundancy_minimisation:
        no_redundancy_path = insert_suffix_before_ext(TRUNCSPEC_OUT_PATH, "_no-redundancy-min")
        write_truncspec(dataset_calls_pre_redundancy, ARC_REF_SEQ_ID, BAC_REF_SEQ_ID, no_redundancy_path)
        print(f"  Wrote: {no_redundancy_path}")
    
    suppress_bootstrap_comments = use_redundancy_minimisation
    
    if use_region_renaming:
        no_renaming_path = insert_suffix_before_ext(TRUNCSPEC_OUT_PATH, "_no-renaming")
        write_truncspec(
            dataset_calls_pre_renaming,
            ARC_REF_SEQ_ID,
            BAC_REF_SEQ_ID,
            no_renaming_path,
            include_bootstrap_comments=not suppress_bootstrap_comments
        )
        print(f"  Wrote: {no_renaming_path}")
    
    write_truncspec(
        dataset_calls_renamed,
        ARC_REF_SEQ_ID,
        BAC_REF_SEQ_ID,
        TRUNCSPEC_OUT_PATH,
        include_bootstrap_comments=not suppress_bootstrap_comments
    )
    print(f"  Wrote: {TRUNCSPEC_OUT_PATH}")
    
    if use_region_renaming:
        renaming_map_path = RESULTS_DIR + "/region_renamings.tsv"
        write_renaming_map(renaming_map, renaming_map_path)
        print(f"  Wrote: {renaming_map_path}")
    print()
    
    # --- 11b. Plot final complete regions ---
    complete_regions = get_complete_regions(dataset_calls_renamed)
    final_region_plots_written = False
    
    if complete_regions:
        print("Plotting final complete regions...")
        plot_final_regions(complete_regions, 'arc', arc_ref_len, FINAL_REGIONS_ARC_PLOT_PATH)
        plot_final_regions(complete_regions, 'bac', bac_ref_len, FINAL_REGIONS_BAC_PLOT_PATH)
        print(f"  Wrote: {FINAL_REGIONS_ARC_PLOT_PATH}")
        print(f"  Wrote: {FINAL_REGIONS_BAC_PLOT_PATH}")
        print()
        final_region_plots_written = True
    else:
        print("No complete regions found; skipping final region plots.")
        print()
    
    # --- 12. Write about file ---
    print("Writing about file...")
    
    config_info = {
        'INFO_OUT_DIR': INFO_OUT_DIR,
        'TRUNCSPEC_OUT_PATH': TRUNCSPEC_OUT_PATH,
        'ARC_REF_SEQ_ID': ARC_REF_SEQ_ID,
        'BAC_REF_SEQ_ID': BAC_REF_SEQ_ID,
        'ARC_HMMER_MAX_EVALUE': ARC_HMMER_MAX_EVALUE,
        'BAC_HMMER_MAX_EVALUE': BAC_HMMER_MAX_EVALUE,
        'MIN_HIT_COVERAGE_FRAC': MIN_HIT_COVERAGE_FRAC,
        'ARC_MIN_HITS_PER_DATASET': ARC_MIN_HITS_PER_DATASET,
        'BAC_MIN_HITS_PER_DATASET': BAC_MIN_HITS_PER_DATASET,
        'OUTLIER_TOL_BP': OUTLIER_TOL_BP,
        'MIN_LEN_BUFFER': MIN_LEN_BUFFER,
        'MAX_LEN_BUFFER': MAX_LEN_BUFFER,
        'USE_CROSS_DOMAIN_BOOTSTRAPPING': USE_CROSS_DOMAIN_BOOTSTRAPPING,
        'CROSS_DOMAIN_BOOTSTRAPPING_MAX_DISTANCE': CROSS_DOMAIN_BOOTSTRAPPING_MAX_DISTANCE,
        'USE_REGION_REDUNDANCY_MINIMISATION': USE_REGION_REDUNDANCY_MINIMISATION,
        'REGION_REDUNDANCY_MINIMISATION_MAX_DISTANCE': REGION_REDUNDANCY_MINIMISATION_MAX_DISTANCE,
        'REGION_REDUNDANCY_MINIMISATION_NAME_JOINER': REGION_REDUNDANCY_MINIMISATION_NAME_JOINER,
        'USE_REGION_RENAMING': USE_REGION_RENAMING,
        'REGION_RENAMING_METHOD': REGION_RENAMING_METHOD
    }
    
    stats_info = {
        'Total ASVs': total_asvs,
        'Total datasets': len(all_datasets),
        'Archaeal hits (after all filtering)': len(arc_gated),
        'Bacterial hits (after all filtering)': len(bac_gated),
        'Archaeal model positions mapped': len(arc_mapping),
        'Bacterial model positions mapped': len(bac_mapping),
        'Datasets with complete regions (pre-redundancy)': complete_count_pre_redundancy,
        'Datasets with incomplete regions (pre-redundancy)': incomplete_count_pre_redundancy,
        'Regions with complete coords (post-redundancy)': complete_count_post_redundancy,
        'Regions with incomplete coords (post-redundancy)': incomplete_count_post_redundancy,
        'Cross-domain bootstrapping applied': use_bootstrapping,
        'Bootstrapped starts': boot_starts,
        'Bootstrapped ends': boot_ends,
        'Region redundancy minimisation applied': use_redundancy_minimisation,
        'Region redundancy merge groups': len(redundancy_merge_groups),
        'Regions removed by redundancy minimisation': len(dataset_calls_pre_redundancy) - len(dataset_calls_redundancy_min),
        'Region renaming applied': use_region_renaming,
        'Regions renamed': len(renaming_map)
    }
    
    about_path = RESULTS_DIR + "/about_processing.txt"
    write_about_file(about_path, config_info, stats_info)
    print(f"  Wrote: {about_path}")
    print()
    
    # --- Done ---
    print("=" * 70)
    print("Step 3 Complete!")
    print("=" * 70)
    print()
    print("Output files written to:")
    print(f"  {TRUNCSPEC_OUT_PATH}")
    if final_region_plots_written:
        print(f"  {FINAL_REGIONS_ARC_PLOT_PATH}")
        print(f"  {FINAL_REGIONS_BAC_PLOT_PATH}")
    print()
    print(f"Intermediate results in {RESULTS_DIR}/:")
    print("  - arc_model_to_ref.tsv")
    print("  - bac_model_to_ref.tsv")
    print("  - asv_hits_arc.tsv")
    print("  - asv_hits_bac.tsv")
    print("  - dataset_region_calls.tsv")
    print("  - dataset_region_qc.tsv")
    print("  - hmmer_filter_stats.tsv")
    if use_bootstrapping:
        print("  - cross_domain_bootstrapping.tsv")
    if use_region_renaming:
        print("  - region_renamings.tsv")
    print("  - about_processing.txt")
    print()
    print(f"Summary: {complete_count_pre_redundancy}/{len(all_datasets)} datasets with complete regions")
    if use_redundancy_minimisation:
        print(f"Final regions: {len(dataset_calls_redundancy_min)} (after redundancy minimisation)")
    if use_region_renaming:
        print(f"Final region names: {len(renaming_map)} regions renamed using '{REGION_RENAMING_METHOD}' method")
    print()
