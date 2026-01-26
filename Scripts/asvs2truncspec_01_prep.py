"""
asvs2truncspec_01_prep.py - Staging/validation script for asvs2truncspec pipeline

This script prepares input data into a clean, predictable structure so that 
later steps (HMMER searches and truncspec generation) don't need to keep
re-discovering metadata.

See asvs2truncspec-plan.md for the overall process and pipeline details.

Processing steps:
1. Create output directory structure
2. Extract reference sequences from the 16S database
3. Load all ASV FASTA files and optional dataset metadata TSVs
4. Validate ASV IDs, dataset IDs, and source filenames
5. Write staged outputs:
   - ref_arc.fna, ref_bac.fna (reference sequences)
   - all_asvs.fna (combined ASVs with standardized headers)
   - dataset_manifest.tsv (ASV-to-dataset mapping)
   - about_staging.txt (summary and configuration snapshot)

Run from the extract16s directory:
    python Scripts/asvs2truncspec_01_prep.py
"""


# Imports
import os
import re
import datetime


# =============================================================================
# CONFIGURATION - Modify these values as needed
# =============================================================================

# Path to the full length 16S GTDB database FASTA file
DB_PATH = "/mnt/secondary/micro16s_dbs/16s_databases/ssu_all_r226.fna"

# Directory containing ASV FASTA files (.fna) and optional dataset TSVs
ASV_DIR = "/home/haig/Repos/extract16s/target_asvs"

# Output directory for info and intermediates
INFO_OUT_DIR = "/home/haig/Repos/extract16s/asvs2truncspec_out/"

# Reference sequence IDs (must exist in DB_PATH)
ARC_REF_SEQ_ID = "RS_GCF_022846175.1~NZ_AP025587.1-#2"
BAC_REF_SEQ_ID = "RS_GCF_030545895.1~NZ_JAUOMX010000042.1"

# Intermediate directory structure (derived from INFO_OUT_DIR)
INTER_DIR = INFO_OUT_DIR + "/intermediates"
STAGED_DIR = INTER_DIR + "/01_staged_inputs"

# Dataset metadata file conventions
DATASET_TSV_SUFFIX = "_datasets.tsv"
DATASET_TSV_ASV_COL = "ASV_ID"
DATASET_TSV_DATASET_COL = "Dataset_ID"
DATASET_TSV_REGION_COL = "Region"  # optional column

# Default region value when no Region metadata is available
UNKNOWN_REGION = "UNK"


# =============================================================================
# VALIDATION HELPERS
# =============================================================================

# Characters that are forbidden in dataset/ASV/source IDs
FORBIDDEN_CHARS_PATTERN = re.compile(r'[|=\s\x00-\x1f\x7f]')


def validate_id(id_str, id_type):
    """Check that an ID doesn't contain forbidden characters.
    
    Forbidden: pipe (|), equals (=), whitespace, non-printing characters.
    Returns error message string if invalid, None if valid.
    """
    if FORBIDDEN_CHARS_PATTERN.search(id_str):
        return f"{id_type} '{id_str}' contains forbidden characters (|, =, whitespace, or non-printing)"
    return None


def get_fasta_base_name(filename):
    """Extract base name from a FASTA filename (remove .fna or .fasta extension)."""
    if filename.endswith(".fna"):
        return filename[:-4]
    elif filename.endswith(".fasta"):
        return filename[:-6]
    return filename


# =============================================================================
# FASTA PARSING
# =============================================================================

def parse_fasta_file(filepath):
    """Parse a FASTA file, yielding (seq_id, sequence) tuples.
    
    seq_id is the first token after '>' (no spaces).
    Sequences may span multiple lines.
    """
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
                # start new sequence (first token after >)
                seq_id = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
    
    # yield last sequence
    if seq_id is not None:
        yield seq_id, ''.join(seq_lines)


def find_ref_sequence(db_path, ref_id):
    """Find and extract a reference sequence from the database.
    
    Searches for the first sequence whose header contains ref_id.
    Returns (full_header, sequence) or (None, None) if not found.
    """
    with open(db_path, 'r') as f:
        current_header = None
        seq_lines = []
        
        for line in f:
            line_stripped = line.strip()
            if not line_stripped:
                continue
            
            if line_stripped.startswith('>'):
                # check if previous sequence was a match
                if current_header is not None and ref_id in current_header:
                    return current_header, ''.join(seq_lines)
                
                # start new sequence
                current_header = line_stripped[1:]  # remove '>'
                seq_lines = []
            else:
                seq_lines.append(line_stripped)
        
        # check last sequence
        if current_header is not None and ref_id in current_header:
            return current_header, ''.join(seq_lines)
    
    return None, None


# =============================================================================
# DATASET METADATA LOADING
# =============================================================================

def load_dataset_tsv(tsv_path):
    """Load a dataset metadata TSV file.
    
    Returns tuple: (asv_to_dataset, asv_to_region)
    - asv_to_dataset: dict mapping ASV_ID -> Dataset_ID
    - asv_to_region: dict mapping ASV_ID -> Region (or None if no Region column)
    """
    asv_to_dataset = {}
    asv_to_region = {}
    has_region_col = False
    
    with open(tsv_path, 'r') as f:
        # parse header to find column indices
        header = f.readline().strip().split('\t')
        try:
            asv_col_idx = header.index(DATASET_TSV_ASV_COL)
            dataset_col_idx = header.index(DATASET_TSV_DATASET_COL)
        except ValueError as e:
            raise ValueError(f"TSV file {tsv_path} missing required columns: {e}")
        
        # check for optional Region column
        region_col_idx = None
        if DATASET_TSV_REGION_COL in header:
            region_col_idx = header.index(DATASET_TSV_REGION_COL)
            has_region_col = True
        
        # parse data rows
        for line_num, line in enumerate(f, start=2):
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) <= max(asv_col_idx, dataset_col_idx):
                raise ValueError(f"TSV file {tsv_path} line {line_num}: not enough columns")
            
            asv_id = parts[asv_col_idx]
            dataset_id = parts[dataset_col_idx]
            asv_to_dataset[asv_id] = dataset_id
            
            # extract region if column exists
            if has_region_col and region_col_idx < len(parts):
                asv_to_region[asv_id] = parts[region_col_idx]
    
    return asv_to_dataset, asv_to_region if has_region_col else None


# =============================================================================
# MAIN SCRIPT
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("asvs2truncspec Step 1: Staging and Validation")
    print("=" * 70)
    print()
    
    # --- 1. Create output directory structure ---
    print("Creating output directory structure...")
    # Only create directories that this script actually populates
    for dir_path in [INFO_OUT_DIR, INTER_DIR, STAGED_DIR]:
        os.makedirs(dir_path, exist_ok=True)
        print(f"  Created: {dir_path}")
    print()
    
    # --- 2. Extract reference sequences from database ---
    print("Extracting reference sequences from database...")
    print(f"  Database: {DB_PATH}")
    
    # find archaeal reference
    arc_header, arc_seq = find_ref_sequence(DB_PATH, ARC_REF_SEQ_ID)
    if arc_seq is None:
        raise ValueError(f"Archaeal reference sequence '{ARC_REF_SEQ_ID}' not found in database")
    print(f"  Found archaeal ref: {ARC_REF_SEQ_ID} ({len(arc_seq)} bp)")
    
    # find bacterial reference
    bac_header, bac_seq = find_ref_sequence(DB_PATH, BAC_REF_SEQ_ID)
    if bac_seq is None:
        raise ValueError(f"Bacterial reference sequence '{BAC_REF_SEQ_ID}' not found in database")
    print(f"  Found bacterial ref: {BAC_REF_SEQ_ID} ({len(bac_seq)} bp)")
    
    # write reference sequences to staged directory
    arc_ref_path = os.path.join(STAGED_DIR, "ref_arc.fna")
    with open(arc_ref_path, 'w') as f:
        f.write(f">{arc_header}\n{arc_seq}\n")
    
    bac_ref_path = os.path.join(STAGED_DIR, "ref_bac.fna")
    with open(bac_ref_path, 'w') as f:
        f.write(f">{bac_header}\n{bac_seq}\n")
    
    print(f"  Wrote: {arc_ref_path}")
    print(f"  Wrote: {bac_ref_path}")
    print()
    
    # --- 3. Discover ASV files and optional dataset TSVs ---
    print(f"Discovering ASV files in {ASV_DIR}...")
    
    # get list of .fna files
    all_files = os.listdir(ASV_DIR)
    fna_files = sorted([f for f in all_files if f.endswith('.fna') or f.endswith('.fasta')])
    
    if not fna_files:
        raise ValueError(f"No .fna or .fasta files found in {ASV_DIR}")
    
    print(f"  Found {len(fna_files)} FASTA files")
    
    # find matching dataset TSV files
    tsv_files = {f for f in all_files if f.endswith(DATASET_TSV_SUFFIX)}
    print(f"  Found {len(tsv_files)} dataset TSV files")
    print()
    
    # --- 4. Load and validate all ASVs ---
    print("Loading and validating ASVs...")
    
    # track all ASVs and their metadata
    # list of (dataset_id, asv_id, source_fna, region, sequence)
    all_asv_records = []
    
    # track ASV IDs globally for duplicate detection
    seen_asv_ids = {}  # asv_id -> source_fna
    
    # track validation errors
    validation_errors = []
    
    # track region metadata availability
    files_with_region_metadata = 0
    
    # process each FASTA file
    for fna_file in fna_files:
        fna_path = os.path.join(ASV_DIR, fna_file)
        fna_base = get_fasta_base_name(fna_file)
        
        print(f"  Processing {fna_file}...")
        
        # validate source filename
        err = validate_id(fna_base, "Source filename")
        if err:
            validation_errors.append(err)
        
        # check for corresponding dataset TSV
        tsv_file = fna_base + DATASET_TSV_SUFFIX
        has_dataset_tsv = tsv_file in tsv_files
        
        # load dataset metadata if available
        asv_to_dataset = {}
        asv_to_region = None
        if has_dataset_tsv:
            tsv_path = os.path.join(ASV_DIR, tsv_file)
            asv_to_dataset, asv_to_region = load_dataset_tsv(tsv_path)
            print(f"    Loaded dataset metadata from {tsv_file}")
            if asv_to_region is not None:
                files_with_region_metadata += 1
                print(f"    Found Region metadata for {len(asv_to_region)} ASVs")
        
        # load ASVs from FASTA
        asv_count = 0
        for asv_id, sequence in parse_fasta_file(fna_path):
            asv_count += 1
            
            # validate ASV ID
            err = validate_id(asv_id, "ASV ID")
            if err:
                validation_errors.append(f"{fna_file}: {err}")
            
            # check for duplicate ASV IDs (global)
            if asv_id in seen_asv_ids:
                validation_errors.append(
                    f"Duplicate ASV ID '{asv_id}' found in {fna_file} "
                    f"(also in {seen_asv_ids[asv_id]})"
                )
            seen_asv_ids[asv_id] = fna_file
            
            # determine dataset ID
            if has_dataset_tsv:
                # must be in the TSV file
                if asv_id not in asv_to_dataset:
                    validation_errors.append(
                        f"{fna_file}: ASV '{asv_id}' present in FASTA but missing from {tsv_file}"
                    )
                    continue
                
                # combine fna base with dataset ID for uniqueness, avoiding duplication
                # if the dataset ID already includes the filename
                dataset_val = asv_to_dataset[asv_id]
                if dataset_val == fna_base or dataset_val.startswith(fna_base + "_"):
                    dataset_id = dataset_val
                else:
                    dataset_id = fna_base + "_" + dataset_val
            else:
                # use FASTA base name as dataset ID
                dataset_id = fna_base
            
            # determine region (use UNK if not available)
            if asv_to_region is not None and asv_id in asv_to_region:
                region = asv_to_region[asv_id]
            else:
                region = UNKNOWN_REGION
            
            # validate dataset ID
            err = validate_id(dataset_id, "Dataset ID")
            if err:
                validation_errors.append(f"{fna_file}: {err}")
            
            # store record (now includes region)
            all_asv_records.append((dataset_id, asv_id, fna_file, region, sequence))
        
        print(f"    Loaded {asv_count} ASVs")
    
    print()
    
    # --- 5. Check for validation errors ---
    if validation_errors:
        print("VALIDATION ERRORS:")
        for err in validation_errors:
            print(f"  - {err}")
        raise ValueError(f"Validation failed with {len(validation_errors)} error(s)")
    
    print(f"Validation passed: {len(all_asv_records)} ASVs from {len(fna_files)} files")
    print()
    
    # --- 6. Write combined ASV FASTA ---
    print("Writing combined ASV FASTA...")
    all_asvs_path = os.path.join(STAGED_DIR, "all_asvs.fna")
    
    with open(all_asvs_path, 'w') as f:
        for dataset_id, asv_id, source_fna, region, sequence in all_asv_records:
            # format: >dataset=XXX|asv=YYY|src=ZZZ.fna
            header = f">dataset={dataset_id}|asv={asv_id}|src={source_fna}"
            f.write(f"{header}\n{sequence}\n")
    
    print(f"  Wrote {len(all_asv_records)} sequences to {all_asvs_path}")
    print()
    
    # --- 7. Write dataset manifest ---
    print("Writing dataset manifest...")
    manifest_path = os.path.join(STAGED_DIR, "dataset_manifest.tsv")
    
    # collect unique datasets and count ASVs per dataset
    dataset_counts = {}
    for dataset_id, asv_id, source_fna, region, sequence in all_asv_records:
        if dataset_id not in dataset_counts:
            dataset_counts[dataset_id] = 0
        dataset_counts[dataset_id] += 1
    
    # count ASVs with known regions
    asvs_with_known_region = sum(1 for _, _, _, r, _ in all_asv_records if r != UNKNOWN_REGION)
    unique_regions = set(r for _, _, _, r, _ in all_asv_records if r != UNKNOWN_REGION)
    
    with open(manifest_path, 'w') as f:
        f.write("dataset_id\tasv_id\tsource_fna\tregion\n")
        for dataset_id, asv_id, source_fna, region, sequence in all_asv_records:
            f.write(f"{dataset_id}\t{asv_id}\t{source_fna}\t{region}\n")
    
    print(f"  Wrote manifest with {len(all_asv_records)} entries")
    print(f"  Unique datasets: {len(dataset_counts)}")
    for ds_id, count in sorted(dataset_counts.items()):
        print(f"    {ds_id}: {count} ASVs")
    print(f"  Region metadata: {asvs_with_known_region}/{len(all_asv_records)} ASVs with known region")
    if unique_regions:
        print(f"  Unique regions: {', '.join(sorted(unique_regions))}")
    print()
    
    # --- 8. Write staging information file ---
    print("Writing staging information file...")
    about_path = os.path.join(STAGED_DIR, "about_staging.txt")
    
    with open(about_path, 'w') as f:
        now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        f.write(f"About Staging (Step 1)\n")
        f.write(f"Generated: {now}\n")
        f.write(f"\n\n")
        
        f.write("Configuration\n")
        f.write("-" * 40 + "\n")
        f.write(f"DB_PATH = {DB_PATH}\n")
        f.write(f"ASV_DIR = {ASV_DIR}\n")
        f.write(f"INFO_OUT_DIR = {INFO_OUT_DIR}\n")
        f.write(f"ARC_REF_SEQ_ID = {ARC_REF_SEQ_ID}\n")
        f.write(f"BAC_REF_SEQ_ID = {BAC_REF_SEQ_ID}\n")
        f.write(f"DATASET_TSV_SUFFIX = {DATASET_TSV_SUFFIX}\n")
        f.write(f"DATASET_TSV_ASV_COL = {DATASET_TSV_ASV_COL}\n")
        f.write(f"DATASET_TSV_DATASET_COL = {DATASET_TSV_DATASET_COL}\n")
        f.write("\n\n")
        
        f.write("Reference Sequences\n")
        f.write("-" * 40 + "\n")
        f.write(f"Archaea Reference:  {ARC_REF_SEQ_ID} (length: {len(arc_seq)} bp)\n")
        f.write(f"Bacteria Reference: {BAC_REF_SEQ_ID} (length: {len(bac_seq)} bp)\n")
        f.write("\n\n")
        
        f.write("Datasets and ASVs\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total ASVs: {len(all_asv_records)}\n")
        f.write(f"Total Datasets: {len(dataset_counts)}\n")
        f.write(f"Total Source FASTA Files: {len(fna_files)}\n")
        f.write(f"Files with Region Metadata: {files_with_region_metadata}\n")
        f.write(f"ASVs with Known Region: {asvs_with_known_region}/{len(all_asv_records)}\n")
        if unique_regions:
            f.write(f"Unique Regions: {', '.join(sorted(unique_regions))}\n")
        f.write("\n")
        
        f.write("ASV Counts per Dataset:\n")
        for ds_id, count in sorted(dataset_counts.items()):
            f.write(f"  {ds_id}: {count}\n")
        f.write("\n")
        
        f.write("Source Files:\n")
        for fna_file in fna_files:
            f.write(f"  {fna_file}\n")
    
    print(f"  Wrote: {about_path}")
    print()
    
    # --- Done ---
    print("=" * 70)
    print("Step 1 Complete!")
    print("=" * 70)
    print()
    print("Staged outputs written to:")
    print(f"  {STAGED_DIR}/")
    print(f"    - ref_arc.fna ({len(arc_seq)} bp)")
    print(f"    - ref_bac.fna ({len(bac_seq)} bp)")
    print(f"    - all_asvs.fna ({len(all_asv_records)} sequences)")
    print(f"    - dataset_manifest.tsv ({len(dataset_counts)} datasets, region column included)")
    print(f"    - about_staging.txt")
    print()
    if asvs_with_known_region > 0:
        print(f"Region metadata: {asvs_with_known_region} ASVs with known region ({len(unique_regions)} unique regions)")
    else:
        print(f"Region metadata: No region data found (all set to '{UNKNOWN_REGION}')")
    print()
    print("Next step: Run asvs2truncspec_02_hmmer.sh")
    print()