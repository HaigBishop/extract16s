"""
create_truncspec.py - Generate .truncspec files from target ASV sequences

This script locates target ASVs inside full archaeal and bacterial reference
sequences, then reduces redundant start/end positions into a representative
set of regions for use with extract16s.

See how-to-create-truncspec.md for the overall process and usage instructions.

Algorithm summary:
1. Load reference sequences (arc/bac) from text files
2. For each ASV in each input FASTA file, find exact or fuzzy matches
3. Cluster hits globally per domain using greedy intersection algorithm
4. Pair arc/bac reduced regions by source file and length
5. Write .truncspec output with computed min/max lengths
6. Generate summary reports and histogram plots

Run from the extract16s directory:
    python Scripts/create_truncspec.py
"""

import os
import re
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt
import edlib


# =============================================================================
# CONFIGURATION - Modify these values as needed
# =============================================================================

# Directory containing target ASV FASTA files (.fna, .fasta)
TARGET_ASV_DIR = "target_asvs"

# Paths to reference sequence files (relative to script working directory)
ARC_REF_PATH = "16s_reference_regions/arc_ref_seq.txt"
BAC_REF_PATH = "16s_reference_regions/bac_ref_seq.txt"

# Output truncspec filename (will be created in current directory)
TRUNCSPEC_OUT_NAME = "new.truncspec"

# Output directory for intermediate files and summaries
OUT_DIR = "create_truncspec_out"

# Tolerance in base pairs for clustering redundant regions
# Hits within this tolerance of each other can be merged into one region
TOL_BP = 15

# Buffer in bp subtracted from minimum region length
MIN_LEN_BUFFER = 50

# Buffer in bp added to maximum region length
MAX_LEN_BUFFER = 50

# Maximum edit distance (mismatches + indels) for fuzzy matching
# Set to 0 to require exact matches and avoid edlib.
MAX_MISMATCHES = 5

# List of input FASTA filenames to skip (exact names, e.g., ["HMC.fna"])
SKIP_INPUT_FILES = [
    # "AGP.fna",
    # "SAR.fna",
    # "HMC.fna", 
    # "CMR-Untrunc_DA.fna",
    # "CMR-Untrunc_SA.fna",
    # "CMR-Untrunc_SP.fna",
    # "CMR-Untrunc_ST.fna",
    # "CMR-Trunc_DA.fna",
    # "CMR-Trunc_SA.fna",
    # "CMR-Trunc_SP.fna",
    # "CMR-Trunc_ST.fna",
]


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def parse_ref_seq(filepath):
    """Parse a reference sequence file, returning (seq_id, sequence).
    
    File format: comment lines starting with #, then >FULL line, then sequence.
    The seq_id is extracted from the first comment line.
    """
    seq_id = None
    sequence = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # first comment line contains the seq_id
            if line.startswith('#') and seq_id is None:
                # extract id: first space-delimited token after #
                match = re.match(r'#\s*(\S+)', line)
                if match:
                    seq_id = match.group(1)
            # sequence header line
            elif line.startswith('>'):
                continue
            # sequence line
            elif seq_id is not None:
                if sequence is None:
                    sequence = line
                else:
                    sequence += line
    
    return seq_id, sequence


def parse_fasta(filepath):
    """Parse a FASTA file, yielding (seq_id, sequence) tuples."""
    seq_id = None
    seq_lines = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # yield previous sequence if exists
                if seq_id is not None:
                    yield seq_id, ''.join(seq_lines)
                # start new sequence
                seq_id = line[1:].split()[0]  # id is first token after >
                seq_lines = []
            else:
                seq_lines.append(line)
    
    # yield last sequence
    if seq_id is not None:
        yield seq_id, ''.join(seq_lines)


def find_asv_in_ref(asv_seq, ref_seq):
    """Find all occurrences of asv_seq in ref_seq.
    
    Returns list of (start, end) tuples (1-based, inclusive).
    """
    hits = []
    start = 0
    asv_len = len(asv_seq)
    
    while True:
        pos = ref_seq.find(asv_seq, start)
        if pos == -1:
            break
        # convert to 1-based inclusive
        hits.append((pos + 1, pos + asv_len))
        start = pos + 1
    
    return hits


def find_asv_in_ref_fuzzy(asv_seq, ref_seq, max_mismatches):
    """Find approximate occurrences of asv_seq in ref_seq using edlib.

    Returns list of (start, end) tuples (1-based, inclusive). edlib returns
    locations for the best alignments within max_mismatches.
    """
    if max_mismatches <= 0:
        return find_asv_in_ref(asv_seq, ref_seq)

    result = edlib.align(asv_seq, ref_seq, mode="HW", task="locations",
                         k=max_mismatches)
    if result["editDistance"] == -1:
        return []

    hits = []
    for start, end in result["locations"]:
        hits.append((start + 1, end + 1))

    return hits


def greedy_cluster(hits, tol_bp):
    """Cluster hits using greedy intersection algorithm.
    
    Args:
        hits: list of (start, end, source_file) tuples
        tol_bp: tolerance in bp for clustering
    
    Returns:
        list of (rep_start, rep_end, source_files_set) tuples
    """
    if not hits:
        return []
    
    # sort by (start, end) for determinism
    sorted_hits = sorted(hits, key=lambda x: (x[0], x[1]))
    
    clusters = []
    
    # initialize first cluster with first hit
    first = sorted_hits[0]
    cur_start_low = first[0] - tol_bp
    cur_start_high = first[0] + tol_bp
    cur_end_low = first[1] - tol_bp
    cur_end_high = first[1] + tol_bp
    cur_sources = {first[2]}
    
    for hit in sorted_hits[1:]:
        start, end, source = hit
        
        # compute proposed intersection
        new_start_low = max(cur_start_low, start - tol_bp)
        new_start_high = min(cur_start_high, start + tol_bp)
        new_end_low = max(cur_end_low, end - tol_bp)
        new_end_high = min(cur_end_high, end + tol_bp)
        
        # check if intersection is valid
        valid = (new_start_low <= new_start_high and
                 new_end_low <= new_end_high and
                 new_start_low <= new_end_high)
        
        if valid:
            # accept hit into current cluster
            cur_start_low = new_start_low
            cur_start_high = new_start_high
            cur_end_low = new_end_low
            cur_end_high = new_end_high
            cur_sources.add(source)
        else:
            # finalize current cluster
            rep_start = cur_start_low
            rep_end = cur_end_low
            clusters.append((rep_start, rep_end, cur_sources))
            
            # start new cluster
            cur_start_low = start - tol_bp
            cur_start_high = start + tol_bp
            cur_end_low = end - tol_bp
            cur_end_high = end + tol_bp
            cur_sources = {source}
    
    # finalize last cluster
    rep_start = cur_start_low
    rep_end = cur_end_low
    clusters.append((rep_start, rep_end, cur_sources))
    
    return clusters


def get_file_stem(filepath):
    """Extract file stem (name without extension) from path."""
    return Path(filepath).stem


# =============================================================================
# MAIN SCRIPT
# =============================================================================

def main():
    # create output directory
    os.makedirs(OUT_DIR, exist_ok=True)

    # load reference sequences
    print("Loading reference sequences...")
    arc_id, arc_seq = parse_ref_seq(ARC_REF_PATH)
    bac_id, bac_seq = parse_ref_seq(BAC_REF_PATH)
    print(f"  Archaeal ref: {arc_id} ({len(arc_seq)} bp)")
    print(f"  Bacterial ref: {bac_id} ({len(bac_seq)} bp)")
    
    # find all target ASV files
    asv_dir = Path(TARGET_ASV_DIR)
    asv_files = sorted([f for f in asv_dir.iterdir()
                        if f.suffix.lower() in ('.fna', '.fasta')])
    print(f"\nFound {len(asv_files)} ASV files in {TARGET_ASV_DIR}/")
    
    # collect all hits per domain
    # hits are (start, end, source_file_stem)
    arc_hits = []
    bac_hits = []
    
    # track per-file statistics
    file_stats = {}  # stem -> {arc_hits, bac_hits, missing, multi, asv_lengths}
    
    # track all ASV lengths per file for plotting
    asv_lengths_by_file = defaultdict(list)
    
    # process each ASV file
    for asv_file in asv_files:
        if asv_file.name in SKIP_INPUT_FILES:
            print(f"\nSkipping {asv_file.name} (listed in SKIP_INPUT_FILES)")
            continue
        stem = get_file_stem(asv_file)
        print(f"\nProcessing {asv_file.name}...")
        
        file_arc_hits = []
        file_bac_hits = []
        missing_asvs = []
        multi_hit_asvs = []
        
        asv_count = 0
        for asv_id, asv_seq in parse_fasta(asv_file):
            asv_count += 1
            asv_len = len(asv_seq)
            asv_lengths_by_file[stem].append(asv_len)
            
            # search in archaeal reference
            if MAX_MISMATCHES == 0:
                arc_matches = find_asv_in_ref(asv_seq, arc_seq)
                bac_matches = find_asv_in_ref(asv_seq, bac_seq)
            else:
                arc_matches = find_asv_in_ref_fuzzy(
                    asv_seq, arc_seq, MAX_MISMATCHES
                )
                bac_matches = find_asv_in_ref_fuzzy(
                    asv_seq, bac_seq, MAX_MISMATCHES
                )
            
            # record hits
            for start, end in arc_matches:
                arc_hits.append((start, end, stem))
                file_arc_hits.append((asv_id, start, end, asv_len))
            
            for start, end in bac_matches:
                bac_hits.append((start, end, stem))
                file_bac_hits.append((asv_id, start, end, asv_len))
            
            # track missing and multi-hit ASVs
            total_matches = len(arc_matches) + len(bac_matches)
            if total_matches == 0:
                missing_asvs.append((asv_id, asv_len))
            if len(arc_matches) > 1 or len(bac_matches) > 1:
                multi_hit_asvs.append((asv_id, len(arc_matches), len(bac_matches)))
        
        # store stats
        file_stats[stem] = {
            'asv_count': asv_count,
            'arc_hits': file_arc_hits,
            'bac_hits': file_bac_hits,
            'missing': missing_asvs,
            'multi': multi_hit_asvs,
            'asv_lengths': asv_lengths_by_file[stem]
        }
        
        print(f"  {asv_count} ASVs, {len(file_arc_hits)} arc hits, "
              f"{len(file_bac_hits)} bac hits, {len(missing_asvs)} missing")
        
        # write per-file hit tables
        with open(os.path.join(OUT_DIR, f"{stem}_hits.tsv"), 'w') as f:
            f.write("asv_id\tdomain\tstart\tend\tasv_len\n")
            for asv_id, start, end, asv_len in file_arc_hits:
                f.write(f"{asv_id}\tarc\t{start}\t{end}\t{asv_len}\n")
            for asv_id, start, end, asv_len in file_bac_hits:
                f.write(f"{asv_id}\tbac\t{start}\t{end}\t{asv_len}\n")
        
        # write missing ASVs
        if missing_asvs:
            with open(os.path.join(OUT_DIR, f"{stem}_missing.tsv"), 'w') as f:
                f.write("asv_id\tasv_len\n")
                for asv_id, asv_len in missing_asvs:
                    f.write(f"{asv_id}\t{asv_len}\n")
    
    # cluster hits globally per domain
    print("\nClustering hits...")
    arc_clusters = greedy_cluster(arc_hits, TOL_BP)
    bac_clusters = greedy_cluster(bac_hits, TOL_BP)
    print(f"  Archaeal: {len(arc_hits)} hits -> {len(arc_clusters)} clusters")
    print(f"  Bacterial: {len(bac_hits)} hits -> {len(bac_clusters)} clusters")
    
    # write cluster mapping
    with open(os.path.join(OUT_DIR, "arc_clusters.tsv"), 'w') as f:
        f.write("cluster_idx\trep_start\trep_end\trep_len\tsource_files\n")
        for i, (start, end, sources) in enumerate(arc_clusters):
            f.write(f"{i}\t{start}\t{end}\t{end-start+1}\t{','.join(sorted(sources))}\n")
    
    with open(os.path.join(OUT_DIR, "bac_clusters.tsv"), 'w') as f:
        f.write("cluster_idx\trep_start\trep_end\trep_len\tsource_files\n")
        for i, (start, end, sources) in enumerate(bac_clusters):
            f.write(f"{i}\t{start}\t{end}\t{end-start+1}\t{','.join(sorted(sources))}\n")
    
    # pair arc/bac regions by source file
    # for each file, collect reduced regions and pair by length
    print("\nPairing regions...")
    
    # build source -> regions mapping
    arc_by_source = defaultdict(list)  # source -> [(start, end, len)]
    bac_by_source = defaultdict(list)
    
    for start, end, sources in arc_clusters:
        region_len = end - start + 1
        for src in sources:
            arc_by_source[src].append((start, end, region_len))
    
    for start, end, sources in bac_clusters:
        region_len = end - start + 1
        for src in sources:
            bac_by_source[src].append((start, end, region_len))
    
    # collect all paired and unpaired entries
    # entry: (name, arc_start, arc_end, bac_start, bac_end, min_len, max_len, is_paired)
    entries = []
    
    # track which regions have been used (for dedup)
    used_regions = set()  # (arc_start, arc_end, bac_start, bac_end) or single domain tuples
    
    all_sources = set(arc_by_source.keys()) | set(bac_by_source.keys())
    
    for src in sorted(all_sources):
        # get regions for this source
        arc_regs = sorted(arc_by_source.get(src, []), key=lambda x: (x[2], x[0], x[1]))
        bac_regs = sorted(bac_by_source.get(src, []), key=lambda x: (x[2], x[0], x[1]))
        
        # pair by index
        n_paired = min(len(arc_regs), len(bac_regs))
        
        for i in range(n_paired):
            arc_start, arc_end, arc_len = arc_regs[i]
            bac_start, bac_end, bac_len = bac_regs[i]
            
            # compute min/max lengths
            min_len = min(arc_len, bac_len) - MIN_LEN_BUFFER
            max_len = max(arc_len, bac_len) + MAX_LEN_BUFFER
            min_len = max(1, min_len)  # ensure positive
            
            region_key = (arc_start, arc_end, bac_start, bac_end)
            entries.append((src, arc_start, arc_end, bac_start, bac_end,
                           min_len, max_len, True, region_key))
        
        # handle unpaired arc regions
        for i in range(n_paired, len(arc_regs)):
            arc_start, arc_end, arc_len = arc_regs[i]
            min_len = max(1, arc_len - MIN_LEN_BUFFER)
            max_len = arc_len + MAX_LEN_BUFFER
            region_key = ('arc', arc_start, arc_end)
            entries.append((src, arc_start, arc_end, None, None,
                           min_len, max_len, False, region_key))
        
        # handle unpaired bac regions
        for i in range(n_paired, len(bac_regs)):
            bac_start, bac_end, bac_len = bac_regs[i]
            min_len = max(1, bac_len - MIN_LEN_BUFFER)
            max_len = bac_len + MAX_LEN_BUFFER
            region_key = ('bac', bac_start, bac_end)
            entries.append((src, None, None, bac_start, bac_end,
                           min_len, max_len, False, region_key))
    
    # deduplicate entries with same region coordinates but different sources
    # group by region_key and merge source names
    region_to_entries = defaultdict(list)
    for entry in entries:
        region_key = entry[8]
        region_to_entries[region_key].append(entry)
    
    # build final deduplicated entries
    final_entries = []
    for region_key, group in region_to_entries.items():
        # merge source names
        sources = sorted(set(e[0] for e in group))
        name_base = '-'.join(sources)
        
        # use first entry for coordinates (all same)
        first = group[0]
        arc_start, arc_end = first[1], first[2]
        bac_start, bac_end = first[3], first[4]
        min_len, max_len = first[5], first[6]
        is_paired = first[7]
        
        final_entries.append((name_base, arc_start, arc_end, bac_start, bac_end,
                             min_len, max_len, is_paired, sources))
    
    # sort for deterministic output: by name, then lengths, then starts
    def entry_sort_key(e):
        name, arc_start, arc_end, bac_start, bac_end, min_len, max_len, is_paired, sources = e
        # use 0 for None values in sorting
        return (name,
                min_len,
                arc_start or 0,
                bac_start or 0)
    
    final_entries.sort(key=entry_sort_key)
    
    # assign global indices
    indexed_entries = []
    for i, entry in enumerate(final_entries):
        name_base = entry[0]
        name = f"{name_base}_{i+1:03d}"
        indexed_entries.append((name,) + entry[1:])
    
    # write truncspec file
    print(f"\nWriting {TRUNCSPEC_OUT_NAME}...")
    with open(TRUNCSPEC_OUT_NAME, 'w') as f:
        # header comments
        f.write("# 16S Region Truncation Specification File (.truncspec)\n")
        f.write("# Generated by create_truncspec.py\n")
        f.write("# This file specifies regions to be extracted from 16S alignments.\n")
        f.write("#\n")
        f.write(f"# Tolerance: {TOL_BP} bp\n")
        f.write(f"# Min length buffer: {MIN_LEN_BUFFER} bp\n")
        f.write(f"# Max length buffer: {MAX_LEN_BUFFER} bp\n")
        f.write(f"# Max mismatches: {MAX_MISMATCHES}\n")
        f.write("\n")
        
        # reference IDs
        f.write(f"ARC_REF_SEQ_ID: {arc_id}\n")
        f.write(f"BAC_REF_SEQ_ID: {bac_id}\n")
        
        # entries
        for entry in indexed_entries:
            name, arc_start, arc_end, bac_start, bac_end, min_len, max_len, is_paired, sources = entry
            
            if is_paired:
                # paired entry
                line = (f"{name}: arc_start={arc_start}, arc_end={arc_end}, "
                       f"bac_start={bac_start}, bac_end={bac_end}, "
                       f"min_len={min_len}, max_len={max_len}")
            elif arc_start is not None:
                # arc-only entry (commented)
                line = (f"#{name}: arc_start={arc_start}, arc_end={arc_end}, "
                       f"bac_start=NA, bac_end=NA, "
                       f"min_len={min_len}, max_len={max_len}")
            else:
                # bac-only entry (commented)
                line = (f"#{name}: arc_start=NA, arc_end=NA, "
                       f"bac_start={bac_start}, bac_end={bac_end}, "
                       f"min_len={min_len}, max_len={max_len}")
            
            f.write(line + "\n")
    
    print(f"  Wrote {len(indexed_entries)} entries "
          f"({sum(1 for e in indexed_entries if e[7])} paired, "
          f"{sum(1 for e in indexed_entries if not e[7])} unpaired)")
    
    # write summary report
    summary_path = os.path.join(OUT_DIR, "summary.txt")
    with open(summary_path, 'w') as f:
        f.write("create_truncspec.py Summary Report\n")
        f.write("=" * 60 + "\n\n")
        
        f.write("Configuration:\n")
        f.write(f"  TARGET_ASV_DIR: {TARGET_ASV_DIR}\n")
        f.write(f"  ARC_REF_PATH: {ARC_REF_PATH}\n")
        f.write(f"  BAC_REF_PATH: {BAC_REF_PATH}\n")
        f.write(f"  TOL_BP: {TOL_BP}\n")
        f.write(f"  MIN_LEN_BUFFER: {MIN_LEN_BUFFER}\n")
        f.write(f"  MAX_LEN_BUFFER: {MAX_LEN_BUFFER}\n")
        f.write(f"  MAX_MISMATCHES: {MAX_MISMATCHES}\n")
        f.write(f"  SKIP_INPUT_FILES: {SKIP_INPUT_FILES}\n\n")
        
        f.write("References:\n")
        f.write(f"  Archaeal: {arc_id} ({len(arc_seq)} bp)\n")
        f.write(f"  Bacterial: {bac_id} ({len(bac_seq)} bp)\n\n")
        
        f.write("Per-file statistics:\n")
        total_asvs = 0
        total_missing = 0
        for stem in sorted(file_stats.keys()):
            stats = file_stats[stem]
            total_asvs += stats['asv_count']
            total_missing += len(stats['missing'])
            f.write(f"  {stem}:\n")
            f.write(f"    ASVs: {stats['asv_count']}\n")
            f.write(f"    Arc hits: {len(stats['arc_hits'])}\n")
            f.write(f"    Bac hits: {len(stats['bac_hits'])}\n")
            f.write(f"    Missing: {len(stats['missing'])}\n")
            f.write(f"    Multi-hit: {len(stats['multi'])}\n")
        
        f.write(f"\nTotals:\n")
        f.write(f"  Total ASVs: {total_asvs}\n")
        f.write(f"  Total missing: {total_missing}\n")
        f.write(f"  Arc clusters: {len(arc_clusters)}\n")
        f.write(f"  Bac clusters: {len(bac_clusters)}\n")
        f.write(f"  Final entries: {len(indexed_entries)}\n")
        f.write(f"    Paired: {sum(1 for e in indexed_entries if e[7])}\n")
        f.write(f"    Unpaired: {sum(1 for e in indexed_entries if not e[7])}\n")
        
        f.write(f"\nFinal regions:\n")
        for entry in indexed_entries:
            name, arc_start, arc_end, bac_start, bac_end, min_len, max_len, is_paired, sources = entry
            status = "paired" if is_paired else "unpaired"
            f.write(f"  {name}: arc=({arc_start},{arc_end}), bac=({bac_start},{bac_end}), "
                   f"len=[{min_len},{max_len}], {status}, sources={sources}\n")
    
    print(f"  Wrote summary to {summary_path}")
    
    # generate plots
    print("\nGenerating plots...")
    
    # collect region lengths per source file for plotting
    region_lengths_by_source = defaultdict(list)
    for entry in indexed_entries:
        name, arc_start, arc_end, bac_start, bac_end, min_len, max_len, is_paired, sources = entry
        # compute actual region length (use average of arc/bac if paired)
        if is_paired:
            arc_len = arc_end - arc_start + 1
            bac_len = bac_end - bac_start + 1
            region_len = (arc_len + bac_len) / 2
        elif arc_start is not None:
            region_len = arc_end - arc_start + 1
        else:
            region_len = bac_end - bac_start + 1
        
        for src in sources:
            region_lengths_by_source[src].append(region_len)
    
    # create histogram for each input file
    for stem in sorted(file_stats.keys()):
        asv_lens = file_stats[stem]['asv_lengths']
        region_lens = region_lengths_by_source.get(stem, [])
        
        if not asv_lens:
            continue
        
        plt.figure(figsize=(10, 6))
        
        # histogram of ASV lengths
        plt.hist(asv_lens, bins=50, alpha=0.7, color='steelblue',
                label=f'ASV lengths (n={len(asv_lens)})')
        
        # overlay region lengths as points
        if region_lens:
            # plot at y=0 with jitter
            y_pos = [0] * len(region_lens)
            plt.scatter(region_lens, y_pos, color='red', s=100, zorder=5,
                        marker='^', label=f'Region lengths (n={len(region_lens)})')
        
        plt.xlabel('Sequence Length (bp)')
        plt.ylabel('Count')
        plt.title(f'{stem}\n({file_stats[stem]["asv_count"]} sequences)')
        plt.legend()
        plt.tight_layout()
        
        plot_path = os.path.join(OUT_DIR, f"{stem}_lengths.png")
        plt.savefig(plot_path, dpi=150)
        plt.close()
    
    print(f"  Saved {len(file_stats)} plots to {OUT_DIR}/")
    
    print("\nDone!")


if __name__ == "__main__":
    main()
