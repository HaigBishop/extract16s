#!/usr/bin/env bash


# extract16s.sh - A script for extracting variable regions from 16S rRNA sequences. Made for Ubuntu.
#  - Uses HMMs to extract regions.
#  - Unless --no_filter_ambiguous flag is given, filter sequences by ambiguous base content (non-ATCG bases)
#  - Unless --no_require_all_regions flag is given, only keep sequences successfully extracted for all regions
#  - If --add_indices flag is given (and --no_require_all_regions is NOT), prepend unique indices (e.g., {1}, {2}) to output FASTA headers.
#  - If using with Micro16S, the <input_fna> given to this script should the output of the filter_database.py script.
#  - If --skip_align flag is given, skips the alignment process and uses pre-existing alignments. 
#       (Assuming they exist at: ./Intermediates/01_align_arc_ssu.a2m and ./Intermediates/01_align_bac_ssu.a2m)
#       (This is good to rerun the truncation process without re-running the alignment process saving much time)
#  - If --rm_intermediates flag is given, removes the Intermediates directory after processing.
#  - If --inter_dir PATH is given, specifies the intermediates directory (default: ./Intermediates)
#  - If --out_dir PATH is given, specifies the output directory (default: ./Output)
#  - If --verbose flag is given, prints verbose output
#  - If --add_indices flag is given, prepends unique indices ({1}, {2}, ...) to output FASTA headers (ignored if --no_require_all_regions is used).



# Input: 16S GTDB database file (e.g. ssu_all_r220_filtered.fna)
#       (You also provide the HMMs and a truncation specification file)


# Output: Directory with results:
#     - FULL_seqs.fasta: Full length sequences
#     - [NAME]_seqs.fasta: Extracted sequences
#     - about_extraction.txt: Summary of processing details


# After running the script you should move the Output directory to somewhere like:
# D:/16S_databases/micro16S_databases/Output/
# (Obviously rename the Output directory to something more meaningful like "m16s_db_007")


# Dependencies:
#  - HMMER suite (hmmalign, esl-reformat)
#  - Unix tools (awk, bc, grep, sed, tr)


# Usage:
#   bash extract16s.sh <input_fna> <bac_hmm> <arc_hmm> <trunc_spec_file> [options]
# 
# Options:
#   --skip_align              Skip alignment process (use existing alignments)
#   --no_filter_ambiguous     Skip filtering by ambiguous base content
#   --no_require_all_regions  Don't require sequences to be present in all regions
#   --rm_intermediates        Remove intermediate files after processing
#   --trunc_padding N         Add N bases of padding to each side of extracted regions (default: 0)
#                            For example, with --trunc_padding 10 and a region specified as
#                            positions 500-700, the actual extraction would be from 490-710
#   --inter_dir PATH          Specify the intermediates directory (default: ./Intermediates)
#   --out_dir PATH            Specify the output directory (default: ./Output)
#   --verbose                 Print verbose output
#   --add_indices             Prepend unique indices ({1}, {2}, ...) to output FASTA headers (ignored if --no_require_all_regions is used).
#
# Example Usage:
#   bash ./Scripts/extract16s.sh \
#      ./InputData/ssu_all_r220_filtered.fna \
#      ./InputData/bac_16s.hmm \
#      ./InputData/arc_16s.hmm \
#      ./InputData/trunc_spec_v4_v3-v4.truncspec \
#      --trunc_padding 15 --rm_intermediates #--skip_align 
#                                                 ^---- OPTIONALLY SKIP ALIGNMENT PROCESS 
#                                                 (if you already have the .a2m alignments)


# Directory structure example:
#  nhmmer_root/
#  │
#  ├── InputData
#  │   ├── arc_16s.hmm                      # HMM database for rRNA models
#  │   ├── bac_16s.hmm                      # HMM database for rRNA models
#  │   ├── ssu_all_r220_filtered.fna        # FASTA file with sequences to analyze
#  │   └── trunc_spec_v4_v3-v4.truncspec    # Truncation specification file
#  │
#  ├── Scripts/
#  │   └── extract16s.sh  # Main script
#  │
#  └── Intermediates/    # (created during script execution, configurable via --inter_dir)
#  │   ├── 00_full_seqs_arc.fna                # Archaea full length sequences
#  │   ├── 00_full_seqs_bac.fna                # Bacteria full length sequences
#  │   ├── 01_align_arc_ssu.sto                # Archaea alignment
#  │   ├── 01_align_bac_ssu.sto                # Bacteria alignment
#  │   ├── 01_align_arc_ssu.a2m                # Archaea alignment (A2M format)
#  │   ├── 01_align_bac_ssu.a2m                # Bacteria alignment (A2M format)
#  │   ├── 02_truncated_arc_{REGION_NAME}.fna        # Archaea truncated sequences
#  │   ├── 02_truncated_bac_{REGION_NAME}.fna        # Bacteria truncated sequences
#  │   ├── 03_joined_truncated_{REGION_NAME}.fna   # joined truncated sequences
#  │   ├── 04_reordered_truncated_{REGION_NAME}.fna  # Reordered truncated sequences
#  │   ├── 05_filtered_truncated_{REGION_NAME}.fna  # Filtered truncated sequences
#  │   └── 06_filtered_full_seqs.fna           # Filtered full sequences
#  │
#  └── Output/           # (created during script execution, configurable via --out_dir)
#      ├── FULL_seqs.fasta                # Full length sequences (same as input FASTA file, but filtered)
#      ├── {REGION_NAME}_seqs.fasta       # Extracted sequences for each region
#      ├── failed_FULL_seqs.fasta         # Sequences that failed filtering
#      └── about_extraction.txt           # Summary of processing details


# More example usage:

# # Ensure all files are in unix format
# find . -type f -name "*.fna" -exec dos2unix {} +
# find . -type f -name "*.hmm" -exec dos2unix {} +
# find . -type f -name "*.sh" -exec dos2unix {} +
# find . -type f -name "*.truncspec" -exec dos2unix {} +

# # Run the script
# bash ./Scripts/extract16s.sh \
#      ./InputData/ssu_all_r220_filtered.fna \
#      ./InputData/bac_16s.hmm \
#      ./InputData/arc_16s.hmm \
#      ./InputData/trunc_spec_v4_v3-v4.truncspec \
#      --trunc_padding 15 --rm_intermediates #--skip_align 
#                                                 ^---- OPTIONALLY SKIP ALIGNMENT PROCESS 
#                                                 (if you already have the .a2m alignments)


# Script Processing Steps
# 1. Split the input fasta file into bacterial and archaeal sequences
# 2. Run hmmalign on the bacterial and archaeal sequences separately
# 3. Convert the alignments to A2M format
# 4. Extract the 16S regions from the A2M alignments based on trunc_padding and the truncation specification file (.truncspec)
#      (the .truncspec file contains reference sequence IDs - each of these strings must be present in a FASTA header of the input file.)
# 5. For each region, write the extracted 16S regions to a new fasta file 
#      (bacteria and archaea merged and in original order with exact same FASTA headers as the input file)
#      e.g. ./Intermediates/ssu_all_r220_filtered_16s_v3.fna
# 6. Filter the sequences and write to the output directory (default: ./Output/)
#      - Filtering by length and ambiguous base content. Then filter for sequences present in all Regions



# Example .truncspec file:
# ```
# ARC_REF_SEQ_ID: GB_GCA_002508705.1~DAYL01000029.1
# BAC_REF_SEQ_ID: GB_GCA_002763075.1~PFDB01000015.1
# V3: arc_start=512, arc_end=775, bac_start=346, bac_end=574, min_len=250, max_len=300
# V3-V4: arc_start=512, arc_end=967, bac_start=346, bac_end=856, min_len=385, max_len=445
# V6: arc_start=1180, arc_end=1432, bac_start=1023, bac_end=1275, min_len=100, max_len=130
# V1-V3: arc_start=45, arc_end=775, bac_start=20, bac_end=574, min_len=450, max_len=580
# ```
# (There is always a ARC_REF_SEQ_ID and a BAC_REF_SEQ_ID. The rest are the regions to be truncated.)
# (This would target 4 regions: V3, V3-V4, V6, and V1-V3. Truncating at positions according to the ref seq IDs)
# (Each region has min_len and max_len parameters to define expected sequence length ranges after truncation)
# If the input FASTA was ./InputData/ssu_all_r220_filtered.fna the output files would be:
#     ./Output/FULL_seqs.fasta
#     ./Output/V3_seqs.fasta
#     ./Output/V3-V4_seqs.fasta
#     ./Output/V6_seqs.fasta
#     ./Output/V1-V3_seqs.fasta
#     ./Output/about_extraction.txt  


# Example input .fna file (sequences shortened for brevity):
# >GB_GCA_000488806.1~KI535272.1 d__Bacteria;p__Bacillota_A;c__Clostridia;o__Peptostreptococcales;f__Anaerovoracaceae;g__Gallibacter;s__Gallibacter brachus [location=411..1920] [ssu_len=1510] [contig_len=1959]
# GGCAGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAACGATGAAGGCCTTTGGGTCGTAAAGTTCTGTTCTAGGTGATGAAAACTGACAGTAACCTAGGAGAAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTACGTAGGTGGCCTT
# >GB_GCA_000492175.2~CP097573.1 d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Lawsonibacter;s__Lawsonibacter sp000492175 [location=1161106..1162683] [ssu_len=1529] [contig_len=3065928]
# GGAGGCAGCAGTGGGGAATATTGGGCAATGGGCGCAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCTTTCGGGTTGTAAACTTCTTTTCTCAGGGACGAAGCAAGTGACGGTACCTGAGGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGG
# >GB_GCA_000492175.2~CP097573.1-#2 d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Lawsonibacter;s__Lawsonibacter sp000492175 [location=2112890..2114418] [ssu_len=1529] [contig_len=3065928]
# GGAGGCAGCAGTGGGGAATATTGGGCAATGGGCGCAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCTTTCGGGTTGTAAACTTCTTTTCTCAGGGACGAAGCAAGTGACGGTACCTGAGGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGG
# >GB_GCA_000493945.1~AWNW01000055.1 d__Bacteria;p__Hydrogenedentota;c__Hydrogenedentia;o__Hydrogenedentiales;f__Hydrogenedentaceae;g__Hydrogenedens;s__Hydrogenedens terephthalicus [location=38027..39574] [ssu_len=1548] [contig_len=80659]
# TAGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGAGGGATGAAGCCCTTCGGGGTGTAAACCTCTTTTCTGGGGGATGAAAAAGATGAGTAGGAAATGACTCATCCTTGACAGTACCCCAGGAATAAGGAACGGCTAACTCCGTGCCAGCAGCCGCGGTAAGACGGAGGTTCCAAGCGTTGTTCGGATTGACTGGGCGTAAAGGGAGCGCAGGCGGTTGAG
# >GB_GCA_000493965.1~AWNV01000024.1 d__Bacteria;p__Acidobacteriota;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__Aminicenans;s__Aminicenans sakinawicola [location=16870..18417] [ssu_len=1548] [contig_len=42608]
# AGCAGTGGGGAATTTTGCGCAATGGGCGAAAGCCTGACGCAGCGACGCCGCGTGGAGGATGAAGGCCTTCGGGTTGTAAACTCCTGTCAGAGGAGAAGAATCCCCGAGTAATCGGGGTTGACGGTATCCTCAAAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCGAACGTTGCTCGGAATTACTGGGCGTAAAGGGTGCGTAGGTGGCTGAGTA
# >GB_GCA_000494145.1~AWOE01000013.1 d__Archaea;p__Thermoproteota;c__Nitrososphaeria_A;o__Caldarchaeales;f__Calditenuaceae;g__Calditenuis;s__Calditenuis sp000494145 [location=10546..12034] [ssu_len=1489] [contig_len=32068]
# CGCCCGTAGCCGGCCCGGTGTGTCCCTCGTTAAATCCACGGGCTTAACCCGTGGGCTGCGGGGGATACTACCGGGCTTGGGGGTGGGAGAGGCGCCCGGTATTCCCGGGGTAGGGGTAAAATCCTCTGATCCCGGGAGGACCATCAGTGGCGAAGGCGGGGCGCCAGAACACGCCCGACGGTGAGGGGCGAAAGCTGGGGGAGCAAACGGGATTAGATACCCCGGTAGTCCCAGCTGTAAACGATGCGGGCTAGCTGTCGGGG
# >GB_GCA_000494185.1~AWOC01000023.1 d__Archaea;p__Thermoproteota;c__Nitrososphaeria_A;o__Caldarchaeales;f__Calditenuaceae;g__AWOC01;s__AWOC01 sp000494185 [location=11503..12993] [ssu_len=1491] [contig_len=14811]



# Example .Intermediates/*.a2m file (sequences shortened for brevity):
# >0000|GB_GCA_000008085.1~AE017199.1 d__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Nanoarchaeales;f__Nanoarchaeaceae;g__Nanoarchaeum;s__Nanoarchaeum equitans [location=432327..433825] [ssu_len=1499] [contig_len=490885]
# c--TCCCGTTGATCCTGCGGGAGGCCACCGCTATCTCCGTCCGGCTAACCCATGGAAGGC
# GAGGGTCCCCGggtaAGGGGGCCCGCCGCACGGCTGAGTAACACGTCGGTAACCTACCCT
# TCAAGCCACCCGAGCTGGGGCCTAGCGAGGCCGTGGGGGGTTcgccCCCCACGGTCGAGC
# TAGGCCCCGGCGAGGGGGGCTAAGTCGACACAAGGTAGCCGTAGGGGAACCTGCGGCTGG
# ATCACCTCC-t
# >0001|GB_GCA_000010605.1~CP005582.1 d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Sulfolobaceae;g__Metallosphaera;s__Metallosphaera sedula [location=1705272..1706706] [ssu_len=1495] [contig_len=2191517]
# a-TTCCGGTTGATCCTGCCGGACCCGATCGCTATAGGGGTAGGGCTAAGCCATGGGAGTC
# GTACGCTCTCGggaAGAGGGCGTGGCGGACGGCTGAGTAACACGTGGCTAACCTGCCCTT
# CACCCGAGTGGAGGGGAAGTGAGGCCTCTTGCCCCTcgggGTGGGAGGTCGAGCTTCTCC
# TCCGCGAGGGGGGAGAAGTCGTAACAAGGTAGCCGTAGGGGAACCTGCGGCTGGATCACC
# TC--
# >0002|GB_GCA_000145985.1~CP002098.1 d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Ignisphaeraceae;g__Ignisphaera;s__Ignisphaera aggregans [location=1279165..1280671] [ssu_len=1507] [contig_len=1875953]
# t--ACCGGTTGATCCTGCCGGACCCGACCGCTATCGGGGTGGGGCTAAGCCATGGAAGTC
# GTACGCCCACCaagtGGTGGGCGTGGCGGACGGCTGAGTAACACGTGGCTAACCTACCCT
# CGGGACGGGGATAGCCCCGGGAAACTGGGGCTAATCCCCGATAGGTGGAGGGGCCTGGAA
# TCGAACCTCTCCTCCGCAAGGGGGGAGAAGTCGTAACAAGGTAGCCGTACCGGAAGGTGC
# GGCTGGATCACCTCC-t
# >0003|GB_GCA_000200715.1~DP000238.1 d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__Cenarchaeum;s__Cenarchaeum symbiosum [location=1200939..1202409] [ssu_len=1471] [contig_len=2045086]
# a--TCCGGTTGATCCTGCCGGACCTGACTGCTATCGGATTGATACTAAGCCATGCGAGTC
# CAGTGAATACGTCCCTGCTCCTTGCACACACCGCCCGTCGTTTCATTGAAGTTTGCTTTT
# AGTGAGGTGACGTCTAAT-TGGCGTTATCGAACTTGAGGTAAGCGACAAGGGAAAAGTCG
# TAACAAGGTGACCGTAGGGGAACCTGCGGTCGGATCACCTCC-t
# >0004|GB_GCA_000204585.1~CM001158.1 d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrososphaerales;f__Nitrosopumilaceae;g__Nitrosarchaeum;s__Nitrosarchaeum limnae [location=573905..575374] [ssu_len=1470] [contig_len=1772718]
# a--TCCGGTTGATCCTGCCGGACCTGACTGCTATCGGAATGATACTAAGCCATGCGAGTC
# ATTGT-AGC----AATACAAGGCAGACGGCTCAGTAACGCGTAGTCAATCTACCCTGTGG
# ACGGGAATAACCTCGGGAAACTGAGAATAATACCCGATAGAACATTATGCCTGGAATGGT
# GCGAGGTGACGTCTGAT-TGGCGTTATCGAACTTGAGGTAAGTGACAAGGGAAAAGTCGT
# AACAAGGTGACCGTAGGGGAACCTGCGGTCGGATCACCTCC-t




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
  local skip_align=false
  
  # Check if --skip_align is among the arguments
  for arg in "$@"; do
    if [ "$arg" = "--skip_align" ]; then
      skip_align=true
      break
    fi
  done
  
  # Check for HMMER tools only if not skipping alignment
  if [ "$skip_align" = false ]; then
    if ! command -v hmmalign >/dev/null 2>&1; then
      missing_deps+=("hmmalign (part of HMMER suite)")
    fi
    
    if ! command -v esl-reformat >/dev/null 2>&1; then
      missing_deps+=("esl-reformat (part of HMMER suite)")
    fi
  fi
  
  # Check for required Unix tools
  for cmd in awk bc grep sed tr; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
      missing_deps+=("$cmd")
    fi
  done
  
  # If any dependencies are missing, print error and exit
  if [ ${#missing_deps[@]} -gt 0 ]; then
    echo "Error: The following required tools are missing:" >&2
    for dep in "${missing_deps[@]}"; do
      echo "  - $dep" >&2
    done
    echo "" >&2
    echo "Please install the missing dependencies before running this script." >&2
    if [ "$skip_align" = false ]; then
      echo "HMMER tools can be obtained from: http://hmmer.org/" >&2
    fi
    exit 1
  fi
}

# Run dependency check
verbose_echo "Checking for required dependencies..."
check_dependencies "$@"




# ===================================================================================================
# Begin Script
# ===================================================================================================
echo "Starting extract16s..."
verbose_echo "Using the following parameters..."
# Start timing the execution
start_time=$(date +%s)





# ===================================================================================================
# Parse Inputs
# ===================================================================================================

# Capture command-line arguments ---------------------------
input_fna="$1"
bac_hmm="$2"
arc_hmm="$3"
trunc_spec_file="$4"
skip_align=false
no_filter_ambiguous=false
no_require_all_regions=false
rm_intermediates=false
trunc_padding=0
verbose=false
out_dir="./Output"
add_indices=false
intermediates_dir="./Intermediates"

i=5  # Start after the four required positional arguments
while [ $i -le $# ]; do
  arg="${!i}"
  
  case "$arg" in
    "--skip_align")
      skip_align=true
      verbose_echo "Skipping alignment process as --skip_align flag is set."
      ;;
    "--no_filter_ambiguous")
      no_filter_ambiguous=true
      verbose_echo "Skipping ambiguous base content filtering as --no_filter_ambiguous flag is set."
      ;;
    "--no_require_all_regions")
      no_require_all_regions=true
      verbose_echo "Skipping requirement for sequences to be present in all regions as --no_require_all_regions flag is set."
      ;;
    "--rm_intermediates")
      rm_intermediates=true
      verbose_echo "Will remove intermediate files after processing as --rm_intermediates flag is set."
      ;;
    "--trunc_padding")
      # Get the next argument as the padding value
      i=$((i+1))  # Move to the next argument instead of using shift
      if [ $i -le $# ]; then  # Make sure there is a next argument
        trunc_padding="${!i}"
        # Validate that it's a non-negative number
        if ! [[ "$trunc_padding" =~ ^[0-9]+$ ]]; then
          echo "Error: trunc_padding must be a non-negative integer"
          exit 1
        fi
        verbose_echo "Using truncation padding of $trunc_padding bases"
      else
        echo "Error: --trunc_padding requires a value"
        exit 1
      fi
      ;;
    "--inter_dir")
      # Get the next argument as the intermediates directory path
      i=$((i+1))
      if [ $i -le $# ]; then
        intermediates_dir="${!i}"
        verbose_echo "Using intermediates directory: $intermediates_dir"
      else
        echo "Error: --inter_dir requires a value"
        exit 1
      fi
      ;;
    "--out_dir")
      # Get the next argument as the output directory path
      i=$((i+1))
      if [ $i -le $# ]; then
        out_dir="${!i}"
        verbose_echo "Using output directory: $out_dir"
      else
        echo "Error: --out_dir requires a value"
        exit 1
      fi
      ;;
    "--verbose")
      verbose=true
      ;;
    "--add_indices")
      add_indices=true
      verbose_echo "Adding unique indices to output FASTA headers as --add_indices flag is set."
      ;;
  esac
  i=$((i+1))
done 

# Check for conflicting options
if [ "$add_indices" = true ] && [ "$no_require_all_regions" = true ]; then
  echo "Warning: --add_indices cannot be used with --no_require_all_regions. Indices will not be added." >&2
  add_indices=false # Prevent adding indices later
fi

verbose_echo "Input parameters:"
verbose_echo "  Input FASTA: $input_fna"
verbose_echo "  Bacteria HMM: $bac_hmm"
verbose_echo "  Archaea HMM: $arc_hmm"
verbose_echo "  Truncation specification file: $trunc_spec_file"
verbose_echo "  Intermediates Directory: $intermediates_dir"
verbose_echo "  Output Directory: $out_dir"


# Validate input files ---------------------------
if [ ! -f "$input_fna" ]; then
  echo "Error: Input FASTA file '$input_fna' not found."
  exit 1
fi

if [ ! -f "$bac_hmm" ]; then
  echo "Error: Bacteria HMM file '$bac_hmm' not found."
  exit 1
fi

if [ ! -f "$arc_hmm" ]; then
  echo "Error: Archaea HMM file '$arc_hmm' not found."
  exit 1
fi

if [ ! -f "$trunc_spec_file" ]; then
  echo "Error: Truncation specification file '$trunc_spec_file' not found."
  exit 1
fi


# Parse the truncation specification file --------------------------- 
verbose_echo "Truncation specifications:"
declare -A region_specs
while IFS= read -r line || [[ -n "$line" ]]; do
  # Skip comments and blank lines
  [[ "$line" =~ ^#.*$ || -z "$line" ]] && continue

  # Clean the line to remove any hidden characters
  line=$(echo "$line" | tr -d '\r')
  
  if [[ "$line" =~ ^ARC_REF_SEQ_ID ]]; then
    arc_ref=$(echo "$line" | cut -d':' -f2 | xargs)
    verbose_echo "  Found Archaea ref sequence ID: $arc_ref"
  elif [[ "$line" =~ ^BAC_REF_SEQ_ID ]]; then
    bac_ref=$(echo "$line" | cut -d':' -f2 | xargs)
    verbose_echo "  Found Bacteria ref sequence ID: $bac_ref"
  elif [[ "$line" =~ : ]]; then
    # More robust parsing of region name and parameters
    region=$(echo "$line" | cut -d':' -f1 | tr -d ' ')
    params=$(echo "$line" | cut -d':' -f2 | xargs)
    
    # Apply padding to start/end positions if trunc_padding > 0
    if [ "$trunc_padding" -gt 0 ]; then
      # Extract current values
      arc_start=$(echo "$params" | grep -oE 'arc_start=[0-9]+' | cut -d'=' -f2)
      arc_end=$(echo "$params" | grep -oE 'arc_end=[0-9]+' | cut -d'=' -f2)
      bac_start=$(echo "$params" | grep -oE 'bac_start=[0-9]+' | cut -d'=' -f2)
      bac_end=$(echo "$params" | grep -oE 'bac_end=[0-9]+' | cut -d'=' -f2)
      
      # Apply padding
      arc_start=$((arc_start - trunc_padding))
      arc_end=$((arc_end + trunc_padding))
      bac_start=$((bac_start - trunc_padding))
      bac_end=$((bac_end + trunc_padding))
      
      # Ensure we don't go below 1
      [ "$arc_start" -lt 1 ] && arc_start=1
      [ "$bac_start" -lt 1 ] && bac_start=1
      
      # Replace values in params string
      params=$(echo "$params" | sed "s/arc_start=[0-9]\+/arc_start=$arc_start/")
      params=$(echo "$params" | sed "s/arc_end=[0-9]\+/arc_end=$arc_end/")
      params=$(echo "$params" | sed "s/bac_start=[0-9]\+/bac_start=$bac_start/")
      params=$(echo "$params" | sed "s/bac_end=[0-9]\+/bac_end=$bac_end/")
    fi
    
    region_specs["$region"]="$params"
    verbose_echo "  Found '$region': $params"
  fi
done < "$trunc_spec_file"


# Validate sequence IDs ---------------------------
if [ -z "$arc_ref" ]; then
  echo "Error: Archaea reference sequence ID not found in truncation specification file."
  exit 1
fi

if [ -z "$bac_ref" ]; then
  echo "Error: Bacteria reference sequence ID not found in truncation specification file."
  exit 1
fi
# Check if reference IDs exist in input FASTA
if ! grep -q "$arc_ref" "$input_fna"; then
  echo "Error: Archaea reference sequence ID '$arc_ref' not found in input FASTA file."
  exit 1
fi

if ! grep -q "$bac_ref" "$input_fna"; then
  echo "Error: Bacteria reference sequence ID '$bac_ref' not found in input FASTA file."
  exit 1
fi





# ===================================================================================================
# Create intermediate directory
# ===================================================================================================
# Validate and create intermediate directory
if [ -e "$intermediates_dir" ] && [ ! -d "$intermediates_dir" ]; then
  echo "Error: Intermediate path '$intermediates_dir' exists but is not a directory."
  exit 1
fi
mkdir -p "$intermediates_dir"
verbose_echo ""
verbose_echo "Ensured intermediate directory exists: $intermediates_dir"




# ===================================================================================================
# Perform alignments
# ===================================================================================================
# If we are NOT skipping alignments
if [ "$skip_align" = false ]; then

  # Split sequences into Bacteria and Archaea ---------------------------
  # Define file paths for split FASTA files
  bac_fna="${intermediates_dir}/00_full_seqs_bac.fna"
  arc_fna="${intermediates_dir}/00_full_seqs_arc.fna"

  verbose_echo ""
  verbose_echo "Splitting sequences into Bacteria and Archaea based on FASTA headers..."
  # Split sequences into Bacteria and Archaea by examining the headers
  # "d__Bacteria" vs "d__Archaea" in the FASTA header lines
  awk 'BEGIN {split_type=""; bac_count=0; arc_count=0}
  {
    if ($0 ~ /^>/) {
      if ($0 ~ /d__Bacteria/) {
        split_type="bac"
        bac_count++
        print $0 >> "'$bac_fna'"
      } else if ($0 ~ /d__Archaea/) {
        split_type="arc"
        arc_count++
        print $0 >> "'$arc_fna'"
      } else {
        # If the header does not match either, choose what you prefer
        split_type=""
      }
    } else {
      if (split_type == "bac") { print $0 >> "'$bac_fna'" }
      if (split_type == "arc") { print $0 >> "'$arc_fna'" }
    }
  }
  END {
    print "Found " bac_count " bacterial sequences and " arc_count " archaeal sequences"
  }' "$input_fna"
  verbose_echo ""

  
  # Align Bacteria ---------------------------
  # Check if files were created and have content
  if [ -s "$bac_fna" ]; then
    echo "Running hmmalign on Bacteria sequences..."

    # Run hmmalign on the Bacteria subset
    hmmalign -o "${intermediates_dir}/01_align_bac_ssu.sto" \
      "$bac_hmm" "$bac_fna"
    echo "Bacteria alignment complete: ${intermediates_dir}/01_align_bac_ssu.sto"
    
    # Convert the alignment to A2M format 
    # (this converts it to DNA and .a2m is much easier to deal with than .sto)
    echo "Converting Bacteria alignment to A2M format..."
    esl-reformat --informat stockholm -d -o "${intermediates_dir}/01_align_bac_ssu.a2m" a2m "${intermediates_dir}/01_align_bac_ssu.sto"
    echo "Bacteria A2M conversion complete: ${intermediates_dir}/01_align_bac_ssu.a2m"
  else
    echo "Warning: No bacterial sequences found or file is empty."
  fi

  # Align Archaea ---------------------------
  if [ -s "$arc_fna" ]; then
    echo "Running hmmalign on Archaea sequences..."

    # Run hmmalign on the Archaea subset
    hmmalign -o "${intermediates_dir}/01_align_arc_ssu.sto" \
      "$arc_hmm" "$arc_fna"
    echo "Archaea alignment complete: ${intermediates_dir}/01_align_arc_ssu.sto"
    
    # Convert the alignment to A2M format 
    # (this converts it to DNA and .a2m is much easier to deal with than .sto)
    echo "Converting Archaea alignment to A2M format..."
    esl-reformat --informat stockholm -d -o "${intermediates_dir}/01_align_arc_ssu.a2m" a2m "${intermediates_dir}/01_align_arc_ssu.sto"
    echo "Archaea A2M conversion complete: ${intermediates_dir}/01_align_arc_ssu.a2m"

  else
    echo "Warning: No archaeal sequences found or file is empty."
  fi
fi




# ===================================================================================================
# Verify A2M alignments
# ===================================================================================================
# Check that the two .a2m files exist
if [ ! -f "${intermediates_dir}/01_align_bac_ssu.a2m" ]; then
  echo "Error: Bacterial alignment file does not exist."
  exit 1
fi
if [ ! -f "${intermediates_dir}/01_align_arc_ssu.a2m" ]; then
  echo "Error: Archaea alignment file does not exist."
  exit 1
fi
# Check that the two .a2m files are not empty
if [ ! -s "${intermediates_dir}/01_align_bac_ssu.a2m" ]; then
  echo "Error: Bacterial alignment file is empty."
  exit 1
fi
if [ ! -s "${intermediates_dir}/01_align_arc_ssu.a2m" ]; then
  echo "Error: Archaea alignment file is empty."
  exit 1
fi




# ===================================================================================================
# Beginning truncation & extraction process
# ===================================================================================================
# Here we loop through each region originating from the truncation specification file:
#   1. Identify the reference sequence and extract that region specified by start and end positions for all sequences for bacteria and archaea separately.
#   2. Merge all the truncated sequences together with the same header and the same order as seen in the original input FASTA file, and write it to file.
#        e.g.   ./InputData/ssu_all_r220_filtered.fna  ->  ./Intermediates/ssu_all_r220_filtered_16s_v3.fna
# This is done for every region (e.g. V3, V3-V4, V6, V1-V3, etc.)
verbose_echo ""
verbose_echo "Beginning extraction of 16S regions using alignments..."





# ===================================================================================================
# Define alignment range helper function
# ===================================================================================================
# Define a helper function to compute the alignment range.

# In the .truncspec file, there are start and end indices which are nucleotide positions in the original reference sequence.
# This function essentially gets the respective alignment indices for the given nucleotide positions.
# This function scans the reference sequence (from the A2M file):
#   - Keeps track of alignment index (only uppercase (A–Z) or the gap symbol ("-"))
#   - Keeps track of the nucleotide index (only uppercase (A–Z) or the lowercase (a–z))
#   - When the nucleotide index equals the requested start or end (from the .truncspec), it records the alignment index.

get_alignment_range() {
  local aln_seq="$1"         # The full aligned reference sequence (a single string)
  local nuc_start_index="$2"  # The start position (from truncspec) in the reference's aligned columns
  local nuc_end_index="$3"    # The end position (from truncspec)
  local pos=0
  local nuc_count=0
  local aln_start=0
  local aln_end=0
  local char
  local aln_seq_len=${#aln_seq}

  for (( i=0; i<aln_seq_len; i++ )); do
    char="${aln_seq:$i:1}"

    # If char is uppercase nucleotide or gap
    if [[ "$char" =~ [A-Z-] ]]; then
      # Count alignment position
      (( pos++ ))
    fi

    # If char is any nucleotide (uppercase or lowercase)
    if [[ "$char" =~ [A-Za-z] ]]; then
      # Count nucleotide
      (( nuc_count++ ))
      # If we have reached the start index, record the alignment index
      if [ "$nuc_count" -eq "$nuc_start_index" ]; then
        aln_start=$pos
      fi
      # If we have reached the end index, record the alignment index and break
      if [ "$nuc_count" -eq "$nuc_end_index" ]; then
        aln_end=$pos
        break
      fi
    fi

  done
  echo "$aln_start $aln_end"
}






# ===================================================================================================
# Region Loop #1
# ===================================================================================================
# For each region (V3, V4, V6, etc.)
for region in "${!region_specs[@]}"; do
  # Skip the ARC_REF_SEQ_ID and BAC_REF_SEQ_ID entries
  if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
    continue
  fi
  verbose_echo ""
  verbose_echo "Processing region: $region"



  # ===================================================================================================
  # Parse .truncspec parameters for the current region
  # ===================================================================================================
  # Parse parameters from the trunc spec string 
  # (e.g. "arc_start=512, arc_end=775, bac_start=346, bac_end=574, min_len=250, max_len=300")
  arc_start=$(echo "${region_specs[$region]}" | grep -oE 'arc_start=[0-9]+' | cut -d'=' -f2)
  arc_end=$(echo "${region_specs[$region]}" | grep -oE 'arc_end=[0-9]+' | cut -d'=' -f2)
  bac_start=$(echo "${region_specs[$region]}" | grep -oE 'bac_start=[0-9]+' | cut -d'=' -f2)
  bac_end=$(echo "${region_specs[$region]}" | grep -oE 'bac_end=[0-9]+' | cut -d'=' -f2)
  min_len=$(echo "${region_specs[$region]}" | grep -oE 'min_len=[0-9]+' | cut -d'=' -f2)
  max_len=$(echo "${region_specs[$region]}" | grep -oE 'max_len=[0-9]+' | cut -d'=' -f2)
  verbose_echo "  Bacteria region positions: $bac_start to $bac_end"
  verbose_echo "  Archaea region positions: $arc_start to $arc_end"
  verbose_echo "  Valid sequence length range: $min_len to $max_len bp"




  # ===================================================================================================
  # Get Bacterial alignment indices
  # ===================================================================================================
  # Here we use Bacterial reference sequence's aligned sequence to compute the alignment indices

  # Find the bacterial reference record (its header should contain the $bac_ref substring)
  ref_bac_record=$(awk -v ref="$bac_ref" 'BEGIN{RS=">"; ORS=""} $0 ~ ref {print ">"$0; exit}' "${intermediates_dir}/01_align_bac_ssu.a2m")
  if [ -z "$ref_bac_record" ]; then
    echo "  Error: Bacterial reference sequence not found in 01_align_bac_ssu.a2m"
    exit 1
  fi
  # Remove the header and newlines to get one continuous sequence string
  ref_bac_seq=$(echo "$ref_bac_record" | sed '1d' | tr -d '\n')
  # Compute the alignment index range for this region based on the Bacterial reference sequence's aligned sequence
  # (use helper function get_alignment_range)
  read bac_aln_start bac_aln_end < <(get_alignment_range "$ref_bac_seq" "$bac_start" "$bac_end")
  verbose_echo "  Bacterial alignment indices: $bac_aln_start to $bac_aln_end"




  # ===================================================================================================
  # Get Archaeal alignment indices
  # ===================================================================================================
  # Here we use Archaeal reference sequence's aligned sequence to compute the alignment indices

  # Find the archaeal reference record (its header should contain the $arc_ref substring)
  ref_arc_record=$(awk -v ref="$arc_ref" 'BEGIN{RS=">"; ORS=""} $0 ~ ref {print ">"$0; exit}' "${intermediates_dir}/01_align_arc_ssu.a2m")
  if [ -z "$ref_arc_record" ]; then
    echo "  Error: Archaeal reference sequence not found in 01_align_arc_ssu.a2m"
    exit 1
  fi
  # Remove the header and newlines to get one continuous sequence string
  ref_arc_seq=$(echo "$ref_arc_record" | sed '1d' | tr -d '\n')
  # Compute the alignment index range for this region based on the Archaeal reference sequence's aligned sequence
  # (use helper function get_alignment_range)
  read arc_aln_start arc_aln_end < <(get_alignment_range "$ref_arc_seq" "$arc_start" "$arc_end")
  verbose_echo "  Archaeal alignment indices: $arc_aln_start to $arc_aln_end"






  # ===================================================================================================
  # Extract truncated regions from the bacterial alignment
  # ===================================================================================================
  # For every record in the Bacterial alignment file (Intermediates/01_align_bac_ssu.a2m) we extract the current region
  # This is done using the Bacterial alignment indices: $bac_aln_start to $bac_aln_end
  # The extracted sequence is written to Intermediates/02_truncated_bac_${region}.fna
  # The alignment indices count uppercase nucleotides and gaps, but NOT lowercase nucleotides.
  # When outputting the sequence, we write all nucleotides (ignore gaps) and convert all lowercase nucleotides to uppercase.
  # We also remove any index prefix (e.g. "112|" or "0002|") from the headers
  # Headers are written as found in the alignment file, without the index prefix.
  
  # Extract truncated sequences to Intermediates/02_truncated_bac_${region}.fna
  truncated_bac_file="${intermediates_dir}/02_truncated_bac_${region}.fna"
  verbose_echo "  Extracting truncated bacterial regions to: $truncated_bac_file"

  awk -v start="$bac_aln_start" -v end="$bac_aln_end" '
  BEGIN {
      RS = ">"
      ORS = ""
  }
  NR > 1 {
      # Split record into header and sequence lines
      pos = index($0, "\n")
      header = substr($0, 1, pos-1)
      seq = substr($0, pos+1)
      gsub(/\n/, "", seq)  # Remove all newlines from sequence
      
      # Remove the index prefix from header
      sub(/^[0-9]+\|/, "", header)
      
      # Extract the region between start and end positions
      # Count only uppercase and gaps for positioning
      current_pos = 0
      extracted = ""
      
      for (i = 1; i <= length(seq); i++) {
          char = substr(seq, i, 1)
          if (char != "." && (char ~ /[A-Za-z-]/)) {
              if (char !~ /[a-z]/) {  # Count uppercase and gaps only
                  current_pos++
              }
              
              # If within our target region, add to extracted sequence
              if (current_pos >= start && current_pos <= end) {
                  if (char != "-") {  # Ignore gaps when extracting
                      # Convert lowercase to uppercase
                      if (char ~ /[a-z]/) {
                          char = toupper(char)
                      }
                      extracted = extracted char
                  }
              }
          }
      }
      
      # Output the header and extracted sequence
      print ">" header "\n" extracted "\n"
  }' "${intermediates_dir}/01_align_bac_ssu.a2m" > "$truncated_bac_file"




  # ===================================================================================================
  # Extract truncated regions from the archaeal alignment
  # ===================================================================================================
  # For every record in the Archaeal alignment file (Intermediates/01_align_arc_ssu.a2m) we extract the current region
  # This is done using the Archaeal alignment indices: $arc_aln_start to $arc_aln_end
  # The extracted sequence is written to Intermediates/02_truncated_arc_${region}.fna
  # The alignment indices count uppercase nucleotides and gaps, but NOT lowercase nucleotides.
  # When outputting the sequence, we write all nucleotides (ignore gaps) and convert all lowercase nucleotides to uppercase.
  # We also remove any index prefix (e.g. "112|" or "0002|") from the headers
  # Headers are written as found in the alignment file, without the index prefix.
  
  # Extract truncated sequences to Intermediates/02_truncated_arc_${region}.fna
  truncated_arc_file="${intermediates_dir}/02_truncated_arc_${region}.fna"
  verbose_echo "  Extracting truncated archaeal regions to: $truncated_arc_file"

  awk -v start="$arc_aln_start" -v end="$arc_aln_end" '
  BEGIN {
      RS = ">"
      ORS = ""
  }
  NR > 1 {
      # Split record into header and sequence lines
      pos = index($0, "\n")
      header = substr($0, 1, pos-1)
      seq = substr($0, pos+1)
      gsub(/\n/, "", seq)  # Remove all newlines from sequence
      
      # Remove the index prefix from header
      sub(/^[0-9]+\|/, "", header)
      
      # Extract the region between start and end positions
      # Count only uppercase and gaps for positioning
      current_pos = 0
      extracted = ""
      
      for (i = 1; i <= length(seq); i++) {
          char = substr(seq, i, 1)
          if (char != "." && (char ~ /[A-Za-z-]/)) {
              if (char !~ /[a-z]/) {  # Count uppercase and gaps only
                  current_pos++
              }
              
              # If within our target region, add to extracted sequence
              if (current_pos >= start && current_pos <= end) {
                  if (char != "-") {  # Ignore gaps when extracting
                      # Convert lowercase to uppercase
                      if (char ~ /[a-z]/) {
                          char = toupper(char)
                      }
                      extracted = extracted char
                  }
              }
          }
      }
      
      # Output the header and extracted sequence
      print ">" header "\n" extracted "\n"
  }' "${intermediates_dir}/01_align_arc_ssu.a2m" > "$truncated_arc_file"





  # ===================================================================================================
  # Merge truncated sequences in original FASTA order
  # ===================================================================================================
  # Now we have two files for this region (e.g. 02_truncated_bac_V4.fna and 02_truncated_arc_V4.fna)
  # We want to put all the truncated sequences in the same order as they appear in the original input FASTA file
  # We do this by joining the two files and then sorting by the original header order
  
  # First, simply concatenate the truncated bacterial and archaeal sequences from 
  # the two files (e.g. 02_truncated_bac_V4.fna and 02_truncated_arc_V4.fna)
  # into a single file (e.g. 03_joined_truncated_V4.fna)
  joined_truncated_file="${intermediates_dir}/03_joined_truncated_${region}.fna"
  cat "$truncated_bac_file" "$truncated_arc_file" > "$joined_truncated_file"

  # Now we will sort the joined truncated sequences by the original header order
  # Output file is Intermediates/04_reordered_truncated_${region}.fna
  merged_output="${intermediates_dir}/04_reordered_truncated_${region}.fna"

  # Merge truncated sequences in original FASTA order.
  # We use the original input FASTA file $input_fna to preserve the header order.
  # First, load all truncated sequences (from both bacteria and archaea) into an associative array keyed by header.
  # Then, for each record in the input FASTA, we extract the header (without the ">" character)
  # and, if a corresponding truncated record exists, print the original header and its truncated sequence.
  awk -v OFS="\n" '
    BEGIN {
      RS=">"; ORS="";
      count = 0;
      match_count = 0;
    }
    # Process the joined truncated sequences file (first input file)
    FNR == NR {
      if ($0 != "") {
        n = split($0, lines, "\n");
        header = lines[1];  # Header (without the leading ">" because RS is ">" )
        # Extract the accession ID (part before first space)
        accession = header;
        if (index(header, " ") > 0) {
          accession = substr(header, 1, index(header, " ")-1);
        }
        seq = "";
        for (i = 2; i <= n; i++) {
          seq = seq lines[i];
        }
        truncated[accession] = seq;  # Save the truncated sequence indexed by ACCESSION ID
        header_map[accession] = header;  # Save the full header for output
        count++;
      }
      next;
    }
    # Process the original FASTA file (second input file) to preserve header order
    {
      if ($0 != "") {
        n = split($0, lines, "\n");
        header = lines[1];  # This header should match one in truncated[]
        # Extract the accession ID (part before first space)
        accession = header;
        if (index(header, " ") > 0) {
          accession = substr(header, 1, index(header, " ")-1);
        }
        if (accession in truncated) {
          # Use the original full header from the input file
          print ">" header "\n" truncated[accession] "\n";
          match_count++;
        }
      }
    }
  ' "$joined_truncated_file" "$input_fna" > "$merged_output"

  verbose_echo "  Region '$region' extraction complete. Output written to: $merged_output"
done




# ===================================================================================================
# Completed extraction for all regions
# ===================================================================================================
# Now for each region we have a file with the truncated sequences of that region in the order of the original input FASTA file
verbose_echo "Completed extraction for all regions."
verbose_echo ""






# ===================================================================================================
# Create and validate output directory
# ===================================================================================================
if [ -e "$out_dir" ] && [ ! -d "$out_dir" ]; then
  echo "Error: Output path '$out_dir' exists but is not a directory."
  exit 1
fi
mkdir -p "$out_dir"
verbose_echo "Ensured output directory exists: $out_dir"
verbose_echo ""





# ===================================================================================================
# Region Loop #2
# ===================================================================================================
# For each region (V3, V4, V6, etc.)
for region in "${!region_specs[@]}"; do
  # Skip the ARC_REF_SEQ_ID and BAC_REF_SEQ_ID entries
  if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
    continue
  fi
  verbose_echo "Filtering truncated sequences for region: $region"



  # ===================================================================================================
  # Filter truncated sequences per region
  # ===================================================================================================
  # For each region, we filter the truncated sequences based on:
  #  - Length must be between (min_len+2*trunc_padding) and (max_len+2*trunc_padding)
  #  - Ambiguous base content must be 0% (unless --no_filter_ambiguous is used)
  # The filtered sequences are written to ./Intermediates/05_filtered_truncated_{REGION_NAME}.fna

  # Extract min and max lengths from the region specification
  min_len=$(echo "${region_specs[$region]}" | grep -oE 'min_len=[0-9]+' | cut -d'=' -f2)
  max_len=$(echo "${region_specs[$region]}" | grep -oE 'max_len=[0-9]+' | cut -d'=' -f2)
  # Compute effective length bounds (each side gets trunc_padding)
  lower_bound=$(( min_len + 2 * trunc_padding ))
  upper_bound=$(( max_len + 2 * trunc_padding ))
  
  # Input and output files
  input_file="${intermediates_dir}/04_reordered_truncated_${region}.fna"
  output_file="${intermediates_dir}/05_filtered_truncated_${region}.fna"
  
  # If ambiguous base filtering is disabled
  if [ "$no_filter_ambiguous" = true ]; then
    # Apply length filtering only
    awk -v lb="$lower_bound" -v ub="$upper_bound" 'BEGIN { RS=">"; ORS="" }
    {
      if($0 != ""){
        n = split($0, lines, "\n")
        header = lines[1]
        seq = ""
        for(i=2; i<=n; i++){
          seq = seq lines[i]
        }
        if(length(seq) >= lb && length(seq) <= ub){
          print ">" header "\n" seq "\n"
        }
      }
    }' "$input_file" > "$output_file"

  # If ambiguous base filtering is enabled
  else
    # Apply both length and ambiguous-base filtering (only A, T, C, G allowed)
    awk -v lb="$lower_bound" -v ub="$upper_bound" 'BEGIN { RS=">"; ORS="" }
    {
      if($0 != ""){
        n = split($0, lines, "\n")
        header = lines[1]
        seq = ""
        for(i=2; i<=n; i++){
          seq = seq lines[i]
        }
        if(length(seq) >= lb && length(seq) <= ub && seq ~ /^[ATCG]+$/){
          print ">" header "\n" seq "\n"
        }
      }
    }' "$input_file" > "$output_file"
  fi
  verbose_echo "  Filtered truncated sequences written to: $output_file"
done







# ===================================================================================================
# Filter full sequences
# ===================================================================================================
# We filter the full length sequences based on:
#  - Ambiguous base content must be 0% (unless --no_filter_ambiguous is used)
# The filtered sequences are written to ./Intermediates/06_filtered_full_seqs.fna

verbose_echo "Filtering full sequences from input FASTA..."

# Output file
full_filtered="${intermediates_dir}/06_filtered_full_seqs.fna"

# If ambiguous base filtering is disabled
if [ "$no_filter_ambiguous" = true ]; then
  # Simply copy the input FASTA
  cp "$input_fna" "$full_filtered"

# If ambiguous base filtering is enabled
else
  # Filter based on ambiguous base content
  awk 'BEGIN { RS=">"; ORS="" }
  {
    if($0 != ""){
      n = split($0, lines, "\n")
      header = lines[1]
      seq = ""
      for(i=2; i<=n; i++){
        seq = seq lines[i]
      }
      # Remove potential problematic characters (like carriage returns) and convert to uppercase
      gsub(/\r/, "", seq)
      seq = toupper(seq)
      # Only print records with sequences of A, T, C, G only (now case-insensitive)
      if(seq ~ /^[ATCG]+$/){
        print ">" header "\n" seq "\n"
      }
    }
  }' "$input_fna" > "$full_filtered"
fi
verbose_echo "  Filtered full sequences written to: $full_filtered"







# ===================================================================================================
# Cross-region filtering
# ===================================================================================================
# Sequences are further filtered ensuring they were successfully truncated and filtered for all regions
# This step is skipped if the --no_require_all_regions flag is given

# If --no_require_all_regions flag is NOT given:
#  - For each sequence, we check if it is present in all regions (including the full sequences)
#    - This is done by looking in files: 
#      - ./Intermediates/05_filtered_truncated_{REGION_NAME}.fna
#      - ./Intermediates/06_filtered_full_seqs.fna
#  - If a sequence is present in all filtered files, it is written to the ./Output/*.fasta files
#  - This means that every ./Output/*.fasta files has the same set of sequences with the same headers in the same order as each other, but different DNA sequences

# If --no_require_all_regions flag is given:
#  - These files are copied (and renamed) to ./Output/
#      - ./Intermediates/05_filtered_truncated_{REGION_NAME}.fna -> ./Output/{REGION}_seqs.fasta
#      - ./Intermediates/06_filtered_full_seqs.fna -> ./Output/FULL_seqs.fasta

# All .fasta output files are written to ./Output/
#  - ./Output/{REGION}_seqs.fasta
#  - ./Output/FULL_seqs.fasta

# If we do not require all regions -----------------------------
if [ "$no_require_all_regions" = true ]; then

  verbose_echo ""
  verbose_echo "Outputting sequences to Output directory..."
  # Copy the filtered full sequences to the Output directory
  cp "$full_filtered" "$out_dir/FULL_seqs.fasta"
  verbose_echo "  Copied filtered full sequences to: $out_dir/FULL_seqs.fasta"

  # Copy the filtered truncated sequences for each region to the Output directory
  for region in "${!region_specs[@]}"; do
    if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
      continue
    fi
    region_filtered="${intermediates_dir}/05_filtered_truncated_${region}.fna"
    cp "$region_filtered" "$out_dir/${region}_seqs.fasta"
    verbose_echo "  Copied filtered truncated sequences for region $region to: $out_dir/${region}_seqs.fasta"
  done



# If we do require all regions -----------------------------
else
  # Filter all the sequences to ensure they are present in all regions and write them in $out_dir/

  # Copy the filtered full sequences to the Output directory
  verbose_echo ""
  verbose_echo "Performing cross-region filtering to retain only sequences present in all filtered files..."


  # Start by extracting accession IDs from the filtered full sequences file
  temp_common="${intermediates_dir}/common_headers.txt"
  grep '^>' "$full_filtered" | sed 's/^>//' | awk '{print $1}' > "$temp_common"
  
  # For each region, intersect its accession IDs with the common accession ID list
  # This creates a progressively smaller list of accession IDs that are present in ALL regions
  for region in "${!region_specs[@]}"; do
    if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
      continue
    fi
    region_file="${intermediates_dir}/05_filtered_truncated_${region}.fna"
    temp_headers="${intermediates_dir}/headers_${region}.txt"
    
    # Extract accession IDs from the current region file
    grep '^>' "$region_file" | sed 's/^>//' | awk '{print $1}' > "$temp_headers"
    
    # Intersect with running list of common accession IDs
    intersect="${intermediates_dir}/common_tmp.txt"
    grep -Fxf "$temp_common" "$temp_headers" > "$intersect"
    mv "$intersect" "$temp_common"
    rm "$temp_headers"
  done
  verbose_echo "  Common headers count: $(wc -l < "$temp_common")"
  
  # Define a helper function to filter a FASTA file by common headers
  # This function takes an input FASTA file and outputs only sequences
  # whose headers are in the common headers list
  filter_by_common_headers() {
    local infile="$1"
    local outfile="$2"
    awk -v common="$temp_common" -v add_indices="$add_indices" 'BEGIN{
      # Read common accession IDs into an array
      while((getline line < common) > 0){
        headers[line]=1
      }
      RS=">"; ORS="";
      idx_count=0 # Initialize index counter
    }
    {
      if($0 != ""){
        n = split($0, lines, "\n")
        header = lines[1]
        # Extract accession ID (part before first space)
        accession = header
        if (index(header, " ") > 0) {
          accession = substr(header, 1, index(header, " ")-1)
        }
        if(accession in headers){
          seq = ""
          for(i=2; i<=n; i++){
            seq = seq lines[i]
          }
          idx_count++ # Increment counter for each sequence written
          # Conditionally prepend index based on add_indices flag
          if (add_indices == "true") {
              print ">{" idx_count "}" header "\n" seq "\n"
          } else {
              print ">" header "\n" seq "\n"
          }
        }
      }
    }' "$infile" > "$outfile"
  }
  
  # Apply the filtering to each file
  # First, filter the full sequences
  filter_by_common_headers "$full_filtered" "$out_dir/FULL_seqs.fasta"
  verbose_echo "  Full sequences after cross-region filtering written to: $out_dir/FULL_seqs.fasta"
  
  # Then for each region, filter its sequences
  for region in "${!region_specs[@]}"; do
    if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
      continue
    fi
    region_file="${intermediates_dir}/05_filtered_truncated_${region}.fna"
    output_region="$out_dir/${region}_seqs.fasta"
    filter_by_common_headers "$region_file" "$output_region"
    verbose_echo "  Filtered truncated sequences for region $region written to: $output_region"
  done

  # Clean up temporary common headers file
  rm "$temp_common"
fi





# ===================================================================================================
# Write failed sequences
# ===================================================================================================
# Write the "failed" sequences that didn't pass filtering
# Specifically, write the full length sequences that are present in the input file but not in the filtered output
verbose_echo "  Writing sequences that failed filtering..."
failed_output="$out_dir/failed_FULL_seqs.fasta"

# Extract accession IDs from the filtered FASTA file into a temporary file, removing any {index} prefix
grep "^>" "$out_dir/FULL_seqs.fasta" | sed 's/^>//' | sed 's/^{.*}//' | awk '{print $1}' > "${intermediates_dir}/passed_headers.txt"

# Compare input against filtered and write sequences not in filtered to failed output
awk -v passed="${intermediates_dir}/passed_headers.txt" 'BEGIN {
  # Load passed accession IDs into array
  while ((getline line < passed) > 0) {
    headers[line] = 1
  }
  RS=">"; ORS=""
}
{
  if ($0 != "") {
    n = split($0, lines, "\n")
    header = lines[1]
    # Extract accession ID (part before first space)
    accession = header
    if (index(header, " ") > 0) {
      accession = substr(header, 1, index(header, " ")-1)
    }
    # If accession NOT in the passed headers, write to failed output
    if (!(accession in headers)) {
      seq = ""
      for (i=2; i<=n; i++) {
        seq = seq lines[i]
      }
      print ">" header "\n" seq "\n"
    }
  }
}' "$input_fna" > "$failed_output"





# ===================================================================================================
# Generate about file
# ===================================================================================================
# Writes a summary file to $out_dir/about_extraction.txt
# Contains:
#  - Input parameters and options used
#  - Reference sequence IDs and region parameters
#  - Sequence counts and filtering statistics
#  - Processing details and completion time
verbose_echo ""
verbose_echo "Generating about_extraction.txt summary file..."
{
  echo "16S rRNA Region Extraction Summary"
  echo "=================================="
  echo ""
  echo "Input Parameters:"
  echo "  Input FASTA file: $input_fna"
  echo "  Bacteria HMM: $bac_hmm"
  echo "  Archaea HMM: $arc_hmm"
  echo "  Truncation specification file: $trunc_spec_file"
  echo "  Intermediates directory: $intermediates_dir"
  echo "  Output directory: $out_dir"
  echo "  Options:"
  echo "    Skip alignment: $skip_align"
  if [ "$no_filter_ambiguous" = true ]; then
    echo "    Filter ambiguous bases: Disabled"
  else
    echo "    Filter ambiguous bases: Enabled"
  fi
  if [ "$no_require_all_regions" = true ]; then
    echo "    Require all regions: Disabled"
    echo "    Add indices: Disabled (incompatible with --no_require_all_regions)"
  else
    echo "    Require all regions: Enabled"
    if [ "$add_indices" = true ]; then
      echo "    Add indices: Enabled"
    else
      echo "    Add indices: Disabled"
    fi
  fi
  echo "    Remove intermediates: $rm_intermediates"
  echo "    Truncation padding: $trunc_padding"
  echo ""
  echo "Reference Sequence IDs:"
  echo "  Archaea: $arc_ref"
  echo "  Bacteria: $bac_ref"
  echo ""
  echo "Region Extraction Parameters:"
  for region in "${!region_specs[@]}"; do
    if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
      continue
    fi
    arc_start=$(echo "${region_specs[$region]}" | grep -oE 'arc_start=[0-9]+' | cut -d'=' -f2)
    arc_end=$(echo "${region_specs[$region]}" | grep -oE 'arc_end=[0-9]+' | cut -d'=' -f2)
    bac_start=$(echo "${region_specs[$region]}" | grep -oE 'bac_start=[0-9]+' | cut -d'=' -f2)
    bac_end=$(echo "${region_specs[$region]}" | grep -oE 'bac_end=[0-9]+' | cut -d'=' -f2)
    min_len=$(echo "${region_specs[$region]}" | grep -oE 'min_len=[0-9]+' | cut -d'=' -f2)
    max_len=$(echo "${region_specs[$region]}" | grep -oE 'max_len=[0-9]+' | cut -d'=' -f2)
    effective_min=$(( min_len + 2 * trunc_padding ))
    effective_max=$(( max_len + 2 * trunc_padding ))
    echo "  Region: $region"
    echo "    Bacteria positions: $bac_start to $bac_end"
    echo "    Archaea positions: $arc_start to $arc_end"
    echo "    Valid sequence lengths (without any padding): $min_len to $max_len bp"
    echo "    Valid sequence length (with padding): $effective_min to $effective_max bp"
  done
  echo ""
  total_seqs=$(grep -c "^>" "$input_fna")
  echo "Total sequences in input FASTA: $total_seqs"
  for region in "${!region_specs[@]}"; do
    if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
      continue
    fi
    region_file="${intermediates_dir}/05_filtered_truncated_${region}.fna"
    region_passed=$(grep -c "^>" "$region_file")
    filtered_out=$(( total_seqs - region_passed ))
    perc_passed=$(awk "BEGIN {printf \"%.2f\", ($region_passed/$total_seqs*100)}")
    perc_filtered=$(awk "BEGIN {printf \"%.2f\", ($filtered_out/$total_seqs*100)}")
    echo "Region $region: $region_passed passed, $filtered_out filtered out ($perc_passed% passed, $perc_filtered% filtered)"
  done
  # Calculate final passed count based on the actual output file
  final_output_full="$out_dir/FULL_seqs.fasta"
  if [ -f "$final_output_full" ]; then
    passed_seqs=$(grep -c "^>" "$final_output_full")
  else
    passed_seqs=0 # Handle case where file might not exist if no seqs pass
  fi
  
  if [ "$no_require_all_regions" = false ]; then
    filtered_all=$(( total_seqs - passed_seqs ))
    perc_passed=$(awk "BEGIN {printf \"%.2f\", ($passed_seqs/$total_seqs*100)}")
    perc_filtered=$(awk "BEGIN {printf \"%.2f\", ($filtered_all/$total_seqs*100)}")
    echo "All regions: $passed_seqs passed, $filtered_all filtered ($perc_passed% passed, $perc_filtered% filtered)"
  fi
  echo ""
  echo "Sequence Length Statistics:"
  calc_length_stats() {
    local fasta_file=$1
    awk 'BEGIN {RS=">"; ORS=""} 
    NR > 1 {
      # Split the record by newline
      n = split($0, lines, "\n")
      # Concatenate lines after the header (lines[1])
      seq = ""
      for (i = 2; i <= n; i++) {
        if (lines[i] != "") { # Avoid adding empty lines
          seq = seq lines[i]
        }
      }
      # Print length if sequence is not empty
      if (seq != "") {
          print length(seq) "\n"
      }
    }' "$fasta_file" | sort -n > "${intermediates_dir}/temp_lengths.txt"
    
    local count=$(wc -l < "${intermediates_dir}/temp_lengths.txt")
    if [ "$count" -gt 0 ]; then
      local min=$(head -n 1 "${intermediates_dir}/temp_lengths.txt")
      local max=$(tail -n 1 "${intermediates_dir}/temp_lengths.txt")
      local median_pos=$(( (count + 1) / 2 ))
      local median=$(sed -n "${median_pos}p" "${intermediates_dir}/temp_lengths.txt")
      local p5_pos=$(echo "($count * 0.05) + 0.5" | bc | awk '{printf "%d", $0}')
      [ $p5_pos -lt 1 ] && p5_pos=1
      local p5=$(sed -n "${p5_pos}p" "${intermediates_dir}/temp_lengths.txt")
      local p95_pos=$(echo "($count * 0.95) + 0.5" | bc | awk '{printf "%d", $0}')
      [ $p95_pos -gt $count ] && p95_pos=$count
      local p95=$(sed -n "${p95_pos}p" "${intermediates_dir}/temp_lengths.txt")
      echo "  Min:$min, 5%:$p5, Median:$median, 95%:$p95, Max:$max"
    else
      echo "  No sequences found"
    fi
    rm "${intermediates_dir}/temp_lengths.txt"
  }
  echo -n "Full sequence lengths: "
  calc_length_stats "$out_dir/FULL_seqs.fasta"
  for region in "${!region_specs[@]}"; do
    if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
      continue
    fi
    echo -n "Region $region lengths: "
    calc_length_stats "$out_dir/${region}_seqs.fasta"
  done
  echo -n "Failed sequence lengths: "
  calc_length_stats "$out_dir/failed_FULL_seqs.fasta"
  echo ""
  echo "Processing completed on: $(date)"
} > "$out_dir/about_extraction.txt"
verbose_echo "about_extraction.txt written to: $out_dir/about_extraction.txt"




# ===================================================================================================
# Clean up intermediate files if requested
# ===================================================================================================
verbose_echo ""
if [ "$rm_intermediates" = true ]; then
  verbose_echo "Cleaning up intermediate files..."
  rm -rf "${intermediates_dir}"
  verbose_echo "Intermediate files removed."
else
  verbose_echo "Intermediate files kept in: ${intermediates_dir}"
  verbose_echo "Use --rm_intermediates flag to remove these files automatically."
fi
verbose_echo ""




# ===================================================================================================
# End script
# ===================================================================================================
# Calculate execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))
minutes=$((execution_time / 60))
seconds=$((execution_time % 60))
# Print done and execution time
echo "Done. All output files have been written to: $out_dir"
echo "Total execution time: ${minutes} minutes and ${seconds} seconds."
