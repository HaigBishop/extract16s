#!/usr/bin/env bash


# asvs2truncspec_02_hmmer.sh - HMMER alignment and search script for asvs2truncspec pipeline
#
# This script is Step 2 of the asvs2truncspec pipeline. It runs HMMER tools to:
#   1. Align reference sequences to their HMMs (hmmalign) - produces Stockholm files
#   2. Search ASV sequences against each HMM (nhmmer) - produces tblout files
#
# The key output from the ASV searches is HMM model coordinates (hmmfrom/hmmto) for each
# ASV hit. Model coordinates are stable across taxa and easy to convert into reference-
# sequence coordinates in Step 3.
#
# See asvs2truncspec-plan.md for the overall process and pipeline details.


# Input files (from Step 1 staged outputs):
#   - ${staged_dir}/ref_arc.fna     : Archaeal reference sequence
#   - ${staged_dir}/ref_bac.fna     : Bacterial reference sequence
#   - ${staged_dir}/all_asvs.fna    : Combined ASV sequences with standardized headers


# Output files (all under ${out_dir}/):
#   - arc_ref.sto       : Stockholm alignment of archaeal ref to arc HMM
#   - bac_ref.sto       : Stockholm alignment of bacterial ref to bac HMM
#   - arc.tblout.gz     : nhmmer tblout results for archaeal HMM (gzipped)
#   - bac.tblout.gz     : nhmmer tblout results for bacterial HMM (gzipped)
#   - arc.hmmer.log     : Log from archaeal nhmmer search
#   - bac.hmmer.log     : Log from bacterial nhmmer search
#   - about_alignment.txt : Summary of run metadata


# Dependencies:
#   - HMMER suite (hmmalign, nhmmer)
#   - gzip
#   - grep


# Usage:
#   bash Scripts/asvs2truncspec_02_hmmer.sh [options]
#
# Required Options:
#   --arc_hmm PATH       Path to archaeal 16S HMM file (e.g., InputData/arc_16s.hmm)
#   --bac_hmm PATH       Path to bacterial 16S HMM file (e.g., InputData/bac_16s.hmm)
#   --staged_dir PATH    Path to Step 1 staged inputs directory
#   --out_dir PATH       Path to output directory for HMMER results
#
# Optional:
#   --threads N          Number of threads for nhmmer (default: 1)
#   --evalue E           E-value threshold for nhmmer reporting (default: 1e-3)
#   --verbose            Print verbose output


# Example usage:
#   bash Scripts/asvs2truncspec_02_hmmer.sh \
#     --arc_hmm InputData/arc_16s.hmm \
#     --bac_hmm InputData/bac_16s.hmm \
#     --staged_dir asvs2truncspec_out/intermediates/01_staged_inputs \
#     --out_dir asvs2truncspec_out/intermediates/02_hmmer \
#     --threads 18 \
#     --verbose




# ===================================================================================================
# Define verbose echo function
# ===================================================================================================
verbose_echo() {
  if [ "$verbose" = true ]; then
    echo "$1"
  fi
}




# ===================================================================================================
# Check for required dependencies
# ===================================================================================================
check_dependencies() {
  local missing_deps=()
  
  # Check for HMMER tools
  if ! command -v hmmalign >/dev/null 2>&1; then
    missing_deps+=("hmmalign (part of HMMER suite)")
  fi
  
  if ! command -v nhmmer >/dev/null 2>&1; then
    missing_deps+=("nhmmer (part of HMMER suite)")
  fi
  
  # Check for gzip
  if ! command -v gzip >/dev/null 2>&1; then
    missing_deps+=("gzip")
  fi
  
  # Check for grep
  if ! command -v grep >/dev/null 2>&1; then
    missing_deps+=("grep")
  fi
  
  # If any dependencies are missing, print error and exit
  if [ ${#missing_deps[@]} -gt 0 ]; then
    echo "Error: The following required tools are missing:" >&2
    for dep in "${missing_deps[@]}"; do
      echo "  - $dep" >&2
    done
    echo "" >&2
    echo "Please install the missing dependencies before running this script." >&2
    echo "HMMER tools can be obtained from: http://hmmer.org/" >&2
    exit 1
  fi
}




# ===================================================================================================
# Parse command-line arguments
# ===================================================================================================

# Initialize defaults
arc_hmm=""
bac_hmm=""
staged_dir=""
out_dir=""
threads=1
evalue="1e-3"
verbose=false

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --arc_hmm)
      arc_hmm="$2"
      shift 2
      ;;
    --bac_hmm)
      bac_hmm="$2"
      shift 2
      ;;
    --staged_dir)
      staged_dir="$2"
      shift 2
      ;;
    --out_dir)
      out_dir="$2"
      shift 2
      ;;
    --threads)
      threads="$2"
      shift 2
      ;;
    --evalue)
      evalue="$2"
      shift 2
      ;;
    --verbose)
      verbose=true
      shift
      ;;
    -h|--help)
      echo "Usage: bash $0 [options]"
      echo ""
      echo "Required Options:"
      echo "  --arc_hmm PATH       Path to archaeal 16S HMM file"
      echo "  --bac_hmm PATH       Path to bacterial 16S HMM file"
      echo "  --staged_dir PATH    Path to Step 1 staged inputs directory"
      echo "  --out_dir PATH       Path to output directory for HMMER results"
      echo ""
      echo "Optional:"
      echo "  --threads N          Number of threads for nhmmer (default: 1)"
      echo "  --evalue E           E-value threshold for nhmmer reporting (default: 1e-3)"
      echo "  --verbose            Print verbose output"
      echo "  -h, --help           Show this help message"
      exit 0
      ;;
    *)
      echo "Error: Unknown option '$1'" >&2
      echo "Use --help for usage information." >&2
      exit 1
      ;;
  esac
done




# ===================================================================================================
# Validate required arguments
# ===================================================================================================
errors=()

if [ -z "$arc_hmm" ]; then
  errors+=("Missing required argument: --arc_hmm")
fi

if [ -z "$bac_hmm" ]; then
  errors+=("Missing required argument: --bac_hmm")
fi

if [ -z "$staged_dir" ]; then
  errors+=("Missing required argument: --staged_dir")
fi

if [ -z "$out_dir" ]; then
  errors+=("Missing required argument: --out_dir")
fi

# Print all errors and exit if any
if [ ${#errors[@]} -gt 0 ]; then
  echo "Error: Missing required arguments:" >&2
  for err in "${errors[@]}"; do
    echo "  - $err" >&2
  done
  echo "" >&2
  echo "Use --help for usage information." >&2
  exit 1
fi




# ===================================================================================================
# Validate input files exist
# ===================================================================================================
echo "========================================================================"
echo "asvs2truncspec Step 2: HMMER Alignment and Search"
echo "========================================================================"
echo ""

verbose_echo "Validating input files..."

# Check HMM files
if [ ! -f "$arc_hmm" ]; then
  echo "Error: Archaeal HMM file not found: $arc_hmm" >&2
  exit 1
fi
verbose_echo "  Found archaeal HMM: $arc_hmm"

if [ ! -f "$bac_hmm" ]; then
  echo "Error: Bacterial HMM file not found: $bac_hmm" >&2
  exit 1
fi
verbose_echo "  Found bacterial HMM: $bac_hmm"

# Check staged input files
ref_arc="${staged_dir}/ref_arc.fna"
ref_bac="${staged_dir}/ref_bac.fna"
all_asvs="${staged_dir}/all_asvs.fna"

if [ ! -f "$ref_arc" ]; then
  echo "Error: Archaeal reference sequence not found: $ref_arc" >&2
  echo "Have you run Step 1 (asvs2truncspec_01_prep.py)?" >&2
  exit 1
fi
verbose_echo "  Found archaeal reference: $ref_arc"

if [ ! -f "$ref_bac" ]; then
  echo "Error: Bacterial reference sequence not found: $ref_bac" >&2
  echo "Have you run Step 1 (asvs2truncspec_01_prep.py)?" >&2
  exit 1
fi
verbose_echo "  Found bacterial reference: $ref_bac"

if [ ! -f "$all_asvs" ]; then
  echo "Error: Combined ASV file not found: $all_asvs" >&2
  echo "Have you run Step 1 (asvs2truncspec_01_prep.py)?" >&2
  exit 1
fi
verbose_echo "  Found combined ASVs: $all_asvs"

verbose_echo ""




# ===================================================================================================
# Check dependencies
# ===================================================================================================
verbose_echo "Checking dependencies..."
check_dependencies
verbose_echo "  All dependencies found"
verbose_echo ""




# ===================================================================================================
# Create output directory
# ===================================================================================================
verbose_echo "Creating output directory..."

if [ -e "$out_dir" ] && [ ! -d "$out_dir" ]; then
  echo "Error: Output path exists but is not a directory: $out_dir" >&2
  exit 1
fi

mkdir -p "$out_dir"
verbose_echo "  Created: $out_dir"
verbose_echo ""




# ===================================================================================================
# Print run parameters
# ===================================================================================================
echo "Run Parameters:"
echo "  Archaeal HMM:     $arc_hmm"
echo "  Bacterial HMM:    $bac_hmm"
echo "  Staged Dir:       $staged_dir"
echo "  Output Dir:       $out_dir"
echo "  Threads:          $threads"
echo "  E-value:          $evalue"
echo ""


# Start timing
start_time=$(date +%s)




# ===================================================================================================
# Step 2a: Align reference sequences to HMMs (hmmalign)
# ===================================================================================================
echo "Aligning reference sequences to HMMs..."

# --- Archaeal reference alignment ---
arc_ref_sto="${out_dir}/arc_ref.sto"
echo "  Running hmmalign for archaeal reference..."
verbose_echo "    Command: hmmalign -o ${arc_ref_sto} ${arc_hmm} ${ref_arc}"

hmmalign -o "${arc_ref_sto}" "${arc_hmm}" "${ref_arc}"

if [ $? -ne 0 ]; then
  echo "Error: hmmalign failed for archaeal reference" >&2
  exit 1
fi
echo "  Created: ${arc_ref_sto}"


# --- Bacterial reference alignment ---
bac_ref_sto="${out_dir}/bac_ref.sto"
echo "  Running hmmalign for bacterial reference..."
verbose_echo "    Command: hmmalign -o ${bac_ref_sto} ${bac_hmm} ${ref_bac}"

hmmalign -o "${bac_ref_sto}" "${bac_hmm}" "${ref_bac}"

if [ $? -ne 0 ]; then
  echo "Error: hmmalign failed for bacterial reference" >&2
  exit 1
fi
echo "  Created: ${bac_ref_sto}"

echo ""




# ===================================================================================================
# Step 2b: Verify RF annotation exists in Stockholm files
# ===================================================================================================
# The RF (reference) annotation line is required for Step 3's model→reference mapping.
# It tells us which alignment columns correspond to HMM match states.

echo "Verifying RF annotation in Stockholm files..."

# Check archaeal alignment
if ! grep -q "^#=GC RF" "${arc_ref_sto}"; then
  echo "Error: No #=GC RF annotation found in ${arc_ref_sto}" >&2
  echo "The HMM may be missing the RF annotation required for coordinate mapping." >&2
  exit 1
fi
verbose_echo "  Found RF annotation in: ${arc_ref_sto}"

# Check bacterial alignment
if ! grep -q "^#=GC RF" "${bac_ref_sto}"; then
  echo "Error: No #=GC RF annotation found in ${bac_ref_sto}" >&2
  echo "The HMM may be missing the RF annotation required for coordinate mapping." >&2
  exit 1
fi
verbose_echo "  Found RF annotation in: ${bac_ref_sto}"

echo "  RF annotations verified successfully"
echo ""




# ===================================================================================================
# Step 2c: Search ASVs against HMMs (nhmmer)
# ===================================================================================================
echo "Searching ASVs against HMMs with nhmmer..."

# Count ASVs for progress info
asv_count=$(grep -c "^>" "${all_asvs}")
echo "  Searching ${asv_count} ASV sequences"


# --- Archaeal HMM search ---
arc_tblout="${out_dir}/arc.tblout"
arc_log="${out_dir}/arc.hmmer.log"

echo "  Running nhmmer with archaeal HMM..."
verbose_echo "    Output: ${arc_tblout}"
verbose_echo "    Log: ${arc_log}"
verbose_echo "    Command: nhmmer --tblout ${arc_tblout} --notextw --noali -E ${evalue} --incE ${evalue} --cpu ${threads} ${arc_hmm} ${all_asvs}"

nhmmer \
  --tblout "${arc_tblout}" \
  --notextw \
  --noali \
  -E "${evalue}" \
  --incE "${evalue}" \
  --cpu "${threads}" \
  "${arc_hmm}" \
  "${all_asvs}" \
  > "${arc_log}" 2>&1

if [ $? -ne 0 ]; then
  echo "Error: nhmmer failed for archaeal HMM search" >&2
  echo "Check log file: ${arc_log}" >&2
  exit 1
fi

# Count hits
arc_hit_count=$(grep -cv "^#" "${arc_tblout}" || echo "0")
echo "  Archaeal search complete: ${arc_hit_count} hits"


# --- Bacterial HMM search ---
bac_tblout="${out_dir}/bac.tblout"
bac_log="${out_dir}/bac.hmmer.log"

echo "  Running nhmmer with bacterial HMM..."
verbose_echo "    Output: ${bac_tblout}"
verbose_echo "    Log: ${bac_log}"
verbose_echo "    Command: nhmmer --tblout ${bac_tblout} --notextw --noali -E ${evalue} --incE ${evalue} --cpu ${threads} ${bac_hmm} ${all_asvs}"

nhmmer \
  --tblout "${bac_tblout}" \
  --notextw \
  --noali \
  -E "${evalue}" \
  --incE "${evalue}" \
  --cpu "${threads}" \
  "${bac_hmm}" \
  "${all_asvs}" \
  > "${bac_log}" 2>&1

if [ $? -ne 0 ]; then
  echo "Error: nhmmer failed for bacterial HMM search" >&2
  echo "Check log file: ${bac_log}" >&2
  exit 1
fi

# Count hits
bac_hit_count=$(grep -cv "^#" "${bac_tblout}" || echo "0")
echo "  Bacterial search complete: ${bac_hit_count} hits"

echo ""




# ===================================================================================================
# Step 2d: Compress tblout files
# ===================================================================================================
echo "Compressing tblout files..."

gzip -f "${arc_tblout}"
echo "  Created: ${arc_tblout}.gz"

gzip -f "${bac_tblout}"
echo "  Created: ${bac_tblout}.gz"

echo ""




# ===================================================================================================
# Step 2e: Write run metadata file
# ===================================================================================================
echo "Writing run metadata..."

about_file="${out_dir}/about_alignment.txt"
end_time=$(date +%s)
elapsed=$((end_time - start_time))

# Get tool versions
hmmalign_version=$(hmmalign -h 2>&1 | head -n 2 | tail -n 1 || echo "unknown")
nhmmer_version=$(nhmmer -h 2>&1 | head -n 2 | tail -n 1 || echo "unknown")

cat > "${about_file}" << EOF
About Alignment (Step 2)
Generated: $(date "+%Y-%m-%d %H:%M:%S")
Elapsed Time: ${elapsed} seconds


Input Files
----------------------------------------
Archaeal HMM:       ${arc_hmm}
Bacterial HMM:      ${bac_hmm}
Archaeal Reference: ${ref_arc}
Bacterial Reference: ${ref_bac}
Combined ASVs:      ${all_asvs}


Output Files
----------------------------------------
Archaeal Ref Alignment: ${arc_ref_sto}
Bacterial Ref Alignment: ${bac_ref_sto}
Archaeal tblout:        ${arc_tblout}.gz
Bacterial tblout:       ${bac_tblout}.gz
Archaeal log:           ${arc_log}
Bacterial log:          ${bac_log}


Parameters
----------------------------------------
Threads:    ${threads}
E-value:    ${evalue}


Tool Versions
----------------------------------------
hmmalign:   ${hmmalign_version}
nhmmer:     ${nhmmer_version}


Search Summary
----------------------------------------
Total ASVs searched: ${asv_count}
Archaeal HMM hits:   ${arc_hit_count}
Bacterial HMM hits:  ${bac_hit_count}


Commands Executed
----------------------------------------
# Archaeal reference alignment
hmmalign -o ${arc_ref_sto} ${arc_hmm} ${ref_arc}

# Bacterial reference alignment
hmmalign -o ${bac_ref_sto} ${bac_hmm} ${ref_bac}

# Archaeal ASV search
nhmmer --tblout ${arc_tblout} --notextw -E ${evalue} --incE ${evalue} --cpu ${threads} ${arc_hmm} ${all_asvs}

# Bacterial ASV search
nhmmer --tblout ${bac_tblout} --notextw -E ${evalue} --incE ${evalue} --cpu ${threads} ${bac_hmm} ${all_asvs}
EOF

echo "  Created: ${about_file}"
echo ""




# ===================================================================================================
# Done
# ===================================================================================================
echo "========================================================================"
echo "Step 2 Complete!"
echo "========================================================================"
echo ""
echo "Output files written to: ${out_dir}/"
echo "  - arc_ref.sto"
echo "  - bac_ref.sto"
echo "  - arc.tblout.gz (${arc_hit_count} hits)"
echo "  - bac.tblout.gz (${bac_hit_count} hits)"
echo "  - arc.hmmer.log"
echo "  - bac.hmmer.log"
echo "  - about_alignment.txt"
echo ""
echo "Elapsed time: ${elapsed} seconds"
echo ""
echo "Next step: Run asvs2truncspec_03_agg.py"
echo ""
