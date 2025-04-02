#!/usr/bin/env bash


# extract_16s_regions_hmmalign.sh - A script for extracting variable regions from 16S rRNA sequences. 
#  - Uses HMMs to extract regions.
#  - Unless --no_filter_ambiguous flag is given, filter sequences by ambiguous base content (non-ATCG bases)
#  - Unless --no_require_all_regions flag is given, only keep sequences that pass through all regions
#  - The <input_fna> given to this script should the output of the filter_database.py script.
#  - If --skip_align flag is given, skips the alignment process and uses pre-existing alignments. 
#       (Assuming they exist at: ./Intermediates/align_arc_ssu.a2m and ./Intermediates/align_bac_ssu.a2m)
#       (This is good to rerun the truncation process without re-running the alignment process saving much time)
#  - If --rm_intermediates flag is given, removes the Intermediates directory after processing.



# Input: 16S GTDB database file (e.g. ssu_all_r220_filtered.fna)
#       (You also provide the HMMs and a truncation specification file)


# Output: Directory with results:
#     - FULL_seqs.fasta: Full length sequences
#     - [NAME]_seqs.fasta: Extracted sequences
#     - about.txt: Summary of processing details


# After running the script you should move the Output directory to somewhere like:
# D:/16S_databases/micro16S_databases/Output/
# (Obviousy rename the Output directory to something more meaningful like "m16s_db_007")


# Dependencies:
#  - HMMER suite (hmmalign, esl-reformat)
#  - Unix tools (awk, bc, grep, sed, tr)


# Usage:
#   bash extract_16s_regions_hmmalign.sh <input_fna> <bac_hmm> <arc_hmm> <trunc_spec_file> [options]
# 
# Options:
#   --skip_align              Skip alignment process (use existing alignments)
#   --no_filter_ambiguous     Skip filtering by ambiguous base content
#   --no_require_all_regions  Don't require sequences to be present in all regions
#   --rm_intermediates        Remove intermediate files after processing
#   --trunc_padding N         Add N bases of padding to each side of extracted regions (default: 0)
#                            For example, with --trunc_padding 10 and a region specified as
#                            positions 500-700, the actual extraction would be from 490-710
#
# Example Usage:
#   bash extract_16s_regions_hmmalign.sh ./InputData/ssu_all_r220_filtered.fna ./InputData/bac_16s.hmm ./InputData/arc_16s.hmm ./InputData/trunc_spec_v3_v3-v4.truncspec





# Directory structure example:

#  nhmmer_root/
#  │
#  ├── InputData
#  │   ├── arc_16s.hmm                      # HMM database for rRNA models
#  │   ├── bac_16s.hmm                      # HMM database for rRNA models
#  │   ├── ssu_all_r220_filtered.fna        # FASTA file with sequences to analyze
#  │   └── trunc_spec_v3_v3-v4.truncspec    # Truncation specification file
#  │
#  ├── Scripts/
#  │   └── extract_16s_regions_hmmalign.sh  # Main script
#  │
#  └── Intermediates/    # (all created during script execution)
#  │   ├── full_seqs_arc.fna                # Archaea full length sequences
#  │   ├── full_seqs_bac.fna                # Bacteria full length sequences
#  │   ├── align_arc_ssu.sto                # Archaea alignment
#  │   ├── align_bac_ssu.sto                # Bacteria alignment
#  │   ├── align_arc_ssu.a2m                # Archaea alignment (A2M format)
#  │   ├── align_bac_ssu.a2m                # Bacteria alignment (A2M format)
#  │   ├── truncated_arc_{REGION_NAME}.fna        # Archaea truncated sequences
#  │   ├── truncated_bac_{REGION_NAME}.fna        # Bacteria truncated sequences
#  │   ├── 01_joined_truncated_{REGION_NAME}.fna   # joined truncated sequences
#  │   ├── 02_reordered_truncated_{REGION_NAME}.fna  # Reordered truncated sequences
#  │   ├── 03_filtered_truncated_{REGION_NAME}.fna  # Filtered truncated sequences
#  │   └── filtered_full_seqs.fna           # Filtered full sequences
#  │
#  └── Output/
#      ├── FULL_seqs.fasta                # Full length sequences (same as input FASTA file, but filtered)
#      ├── {REGION_NAME}_seqs.fasta       # Extracted sequences for each region
#      └── about.txt                      # Summary of processing details




# More example usage:

# # Ensure all files are in unix format
# find . -type f -name "*.fna" -exec dos2unix {} +
# find . -type f -name "*.hmm" -exec dos2unix {} +
# find . -type f -name "*.sh" -exec dos2unix {} +
# find . -type f -name "*.truncspec" -exec dos2unix {} +

# # Run the script
# bash ./Scripts/extract_16s_regions_hmmalign.sh \
#      ./InputData/ssu_all_r220_filtered.fna \
#      ./InputData/bac_16s.hmm \
#      ./InputData/arc_16s.hmm \
#      ./InputData/trunc_spec_v3_v3-v4.truncspec \
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
# 6. Filter the sequences and write to ./Output/
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
#     ./Output/about.txt  


# Example input .fna file (sequences shortened for brevity):
# >GB_GCA_000488855.1~KI535272.1 d__Bacteria;p__Bacillota_A;c__Clostridia;o__Peptostreptococcales;f__Anaerovoracaceae;g__Gallibacter;s__Gallibacter brachus [location=411..1920] [ssu_len=1510] [contig_len=1959]
# GGCAGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAACGATGAAGGCCTTTGGGTCGTAAAGTTCTGTTCTAGGTGATGAAAACTGACAGTAACCTAGGAGAAAGCCCCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGGGCAAGCGTTATCCGGAATTATTGGGCGTAAAGGGTACGTAGGTGGCCTT
# >GB_GCA_000492175.2~CP097573.1 d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Lawsonibacter;s__Lawsonibacter sp000492175 [location=1161155..1162683] [ssu_len=1529] [contig_len=3665928]
# GGAGGCAGCAGTGGGGAATATTGGGCAATGGGCGCAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCTTTCGGGTTGTAAACTTCTTTTCTCAGGGACGAAGCAAGTGACGGTACCTGAGGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGG
# >GB_GCA_000492175.2~CP097573.1-#2 d__Bacteria;p__Bacillota_A;c__Clostridia;o__Oscillospirales;f__Oscillospiraceae;g__Lawsonibacter;s__Lawsonibacter sp000492175 [location=2112890..2114418] [ssu_len=1529] [contig_len=3665928]
# GGAGGCAGCAGTGGGGAATATTGGGCAATGGGCGCAAGCCTGACCCAGCAACGCCGCGTGAAGGAAGAAGGCTTTCGGGTTGTAAACTTCTTTTCTCAGGGACGAAGCAAGTGACGGTACCTGAGGAATAAGCCACGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGTGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGCGTGTAGGCGGG
# >GB_GCA_000493945.1~AWNW01000006.1 d__Bacteria;p__Hydrogenedentota;c__Hydrogenedentia;o__Hydrogenedentiales;f__Hydrogenedentaceae;g__Hydrogenedens;s__Hydrogenedens terephthalicus [location=38027..39574] [ssu_len=1548] [contig_len=86659]
# TAGGGAATATTGGACAATGGGCGAAAGCCTGATCCAGCCATGCCGCGTGAGGGATGAAGCCCTTCGGGGTGTAAACCTCTTTTCTGGGGGATGAAAAAGATGAGTAGGAAATGACTCATCCTTGACAGTACCCCAGGAATAAGGAACGGCTAACTCCGTGCCAGCAGCCGCGGTAAGACGGAGGTTCCAAGCGTTGTTCGGATTGACTGGGCGTAAAGGGAGCGCAGGCGGTTGAG
# >GB_GCA_000493965.1~AWNV01000024.1 d__Bacteria;p__Acidobacteriota;c__Aminicenantia;o__Aminicenantales;f__Aminicenantaceae;g__Aminicenans;s__Aminicenans sakinawicola [location=16870..18417] [ssu_len=1548] [contig_len=42608]
# AGCAGTGGGGAATTTTGCGCAATGGGCGAAAGCCTGACGCAGCGACGCCGCGTGGAGGATGAAGGCCTTCGGGTTGTAAACTCCTGTCAGAGGAGAAGAATCCCCGAGTAATCGGGGTTGACGGTATCCTCAAAGGAAGCCCCGGCTAACTACGTGCCAGCAGCCGCGGTAATACGTAGGGGGCGAACGTTGCTCGGAATTACTGGGCGTAAAGGGTGCGTAGGTGGCTGAGTA
# >GB_GCA_000494145.1~AWOE01000013.1 d__Archaea;p__Thermoproteota;c__Nitrososphaeria_A;o__Caldarchaeales;f__Calditenuaceae;g__Calditenuis;s__Calditenuis sp000494145 [location=10546..12034] [ssu_len=1489] [contig_len=32668]
# CGCCCGTAGCCGGCCCGGTGTGTCCCTCGTTAAATCCACGGGCTTAACCCGTGGGCTGCGGGGGATACTACCGGGCTTGGGGGTGGGAGAGGCGCCCGGTATTCCCGGGGTAGGGGTAAAATCCTCTGATCCCGGGAGGACCATCAGTGGCGAAGGCGGGGCGCCAGAACACGCCCGACGGTGAGGGGCGAAAGCTGGGGGAGCAAACGGGATTAGATACCCCGGTAGTCCCAGCTGTAAACGATGCGGGCTAGCTGTCGGGG
# >GB_GCA_000494185.1~AWOC01000023.1 d__Archaea;p__Thermoproteota;c__Nitrososphaeria_A;o__Caldarchaeales;f__Calditenuaceae;g__AWOC01;s__AWOC01 sp000494185 [location=11503..12993] [ssu_len=1491] [contig_len=14811]



# Example .Intermediates/*.a2m file (sequences shortened for brevity):
# >0000|GB_GCA_000008085.1~AE017199.1 d__Archaea;p__Nanoarchaeota;c__Nanoarchaeia;o__Nanoarchaeales;f__Nanoarchaeaceae;g__Nanoarchaeum;s__Nanoarchaeum equitans [location=432327..433825] [ssu_len=1499] [contig_len=490885]
# c--TCCCGTTGATCCTGCGGGAGGCCACCGCTATCTCCGTCCGGCTAACCCATGGAAGGC
# GAGGGTCCCCGggtaAGGGGGCCCGCCGCACGGCTGAGTAACACGTCGGTAACCTACCCT
# TCAAGCCACCCGAGCTGGGGCCTAGCGAGGCCGTGGGGGGTTcgccCCCCACGGTCGAGC
# TAGGCCCCGGCGAGGGGGGCTAAGTCGACACAAGGTAGCCGTAGGGGAACCTGCGGCTGG
# ATCACCTCC-t
# >0001|GB_GCA_000016605.1~CP000682.1 d__Archaea;p__Thermoproteota;c__Thermoprotei_A;o__Sulfolobales;f__Sulfolobaceae;g__Metallosphaera;s__Metallosphaera sedula [location=1705272..1706766] [ssu_len=1495] [contig_len=2191517]
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




# Check for required dependencies
check_dependencies() {
  local missing_deps=()
  
  # Check for HMMER tools
  if ! command -v hmmalign >/dev/null 2>&1; then
    missing_deps+=("hmmalign (part of HMMER suite)")
  fi
  
  if ! command -v esl-reformat >/dev/null 2>&1; then
    missing_deps+=("esl-reformat (part of HMMER suite)")
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
    echo "HMMER tools can be obtained from: http://hmmer.org/" >&2
    exit 1
  fi
}

# Run dependency check
echo "Checking for required dependencies..."
check_dependencies



# Start timing the execution
start_time=$(date +%s)

echo "Starting extract_16s_regions_hmmalign.sh script..."

# Capture command-line arguments
input_fna="$1"
bac_hmm="$2"
arc_hmm="$3"
trunc_spec_file="$4"
skip_align=false
no_filter_ambiguous=false
no_require_all_regions=false
rm_intermediates=false
trunc_padding=0

# Process arguments
i=5  # Start after the four required positional arguments
while [ $i -le $# ]; do
  arg="${!i}"
  
  case "$arg" in
    "--skip_align")
      skip_align=true
      echo "Skipping alignment process as --skip_align flag is set."
      ;;
    "--no_filter_ambiguous")
      no_filter_ambiguous=true
      echo "Skipping ambiguous base content filtering as --no_filter_ambiguous flag is set."
      ;;
    "--no_require_all_regions")
      no_require_all_regions=true
      echo "Skipping requirement for sequences to be present in all regions as --no_require_all_regions flag is set."
      ;;
    "--rm_intermediates")
      rm_intermediates=true
      echo "Will remove intermediate files after processing as --rm_intermediates flag is set."
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
        echo "Using truncation padding of $trunc_padding bases"
      else
        echo "Error: --trunc_padding requires a value"
        exit 1
      fi
      ;;
  esac
  i=$((i+1))
done 

echo "Input parameters:"
echo "  Input FASTA: $input_fna"
echo "  Bacteria HMM: $bac_hmm"
echo "  Archaea HMM: $arc_hmm"
echo "  Truncation specification file: $trunc_spec_file"

# After parsing arguments
# Validate input files
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

# Parse the truncation specification file into variables
echo "Truncation specifications:"
declare -A region_specs
while IFS= read -r line || [[ -n "$line" ]]; do
  # Skip comments and blank lines
  [[ "$line" =~ ^#.*$ || -z "$line" ]] && continue

  # Clean the line to remove any hidden characters
  line=$(echo "$line" | tr -d '\r')
  
  if [[ "$line" =~ ^ARC_REF_SEQ_ID ]]; then
    arc_ref=$(echo "$line" | cut -d':' -f2 | xargs)
    echo "  Found Archaea ref sequence ID: $arc_ref"
  elif [[ "$line" =~ ^BAC_REF_SEQ_ID ]]; then
    bac_ref=$(echo "$line" | cut -d':' -f2 | xargs)
    echo "  Found Bacteria ref sequence ID: $bac_ref"
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
    echo "  Found '$region': $params"
  fi
done < "$trunc_spec_file"

# After parsing the truncspec file
if [ -z "$arc_ref" ]; then
  echo "Error: Archaea reference sequence ID not found in truncation specification file."
  exit 1
fi

if [ -z "$bac_ref" ]; then
  echo "Error: Bacteria reference sequence ID not found in truncation specification file."
  exit 1
fi

# Define an output directory for intermediate files
intermediates_dir="./Intermediates"
mkdir -p "$intermediates_dir"
echo "Created intermediates directory: $intermediates_dir"



# ------------------------------
# Alignment Process
# ------------------------------
if [ "$skip_align" = false ]; then
  # Define file paths for split FASTA files
  bac_fna="${intermediates_dir}/full_seqs_bac.fna"
  arc_fna="${intermediates_dir}/full_seqs_arc.fna"

  echo "Splitting sequences into Bacteria and Archaea based on FASTA headers..."
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

  # Check if files were created and have content
  if [ -s "$bac_fna" ]; then
    echo "Running hmmalign on Bacteria sequences..."
    # Run hmmalign on the Bacteria subset
    hmmalign -o "${intermediates_dir}/align_bac_ssu.sto" \
      "$bac_hmm" "$bac_fna"
    echo "Bacteria alignment complete: ${intermediates_dir}/align_bac_ssu.sto"
    
    # Convert the alignment to A2M format (DNA and much easier to deal with that .sto)
    echo "Converting Bacteria alignment to A2M format..."
    esl-reformat --informat stockholm -d -o "${intermediates_dir}/align_bac_ssu.a2m" a2m "${intermediates_dir}/align_bac_ssu.sto"
    echo "Bacteria A2M conversion complete: ${intermediates_dir}/align_bac_ssu.a2m"
  else
    echo "Warning: No bacterial sequences found or file is empty."
  fi

  if [ -s "$arc_fna" ]; then
    echo "Running hmmalign on Archaea sequences..."
    # Run hmmalign on the Archaea subset
    hmmalign -o "${intermediates_dir}/align_arc_ssu.sto" \
      "$arc_hmm" "$arc_fna"
    echo "Archaea alignment complete: ${intermediates_dir}/align_arc_ssu.sto"
    
    # Convert the alignment to A2M format (DNA and much easier to deal with that .sto)
    echo "Converting Archaea alignment to A2M format..."
    esl-reformat --informat stockholm -d -o "${intermediates_dir}/align_arc_ssu.a2m" a2m "${intermediates_dir}/align_arc_ssu.sto"
    echo "Archaea A2M conversion complete: ${intermediates_dir}/align_arc_ssu.a2m"

  else
    echo "Warning: No archaeal sequences found or file is empty."
  fi


  # Check that the two .a2m files exist
  if [ ! -f "${intermediates_dir}/align_bac_ssu.a2m" ]; then
    echo "Error: Bacterial alignment file does not exist."
    exit 1
  fi
  if [ ! -f "${intermediates_dir}/align_arc_ssu.a2m" ]; then
    echo "Error: Archaea alignment file does not exist."
    exit 1
  fi
  # Check that the two .a2m files are not empty
  if [ ! -s "${intermediates_dir}/align_bac_ssu.a2m" ]; then
    echo "Error: Bacterial alignment file is empty."
    exit 1
  fi
  if [ ! -s "${intermediates_dir}/align_arc_ssu.a2m" ]; then
    echo "Error: Archaea alignment file is empty."
    exit 1
  fi
fi




# ------------------------------
# Truncation & Extraction Process
# ------------------------------
# We loop through each region originating from the truncation specification file:
#   1. Identify the reference sequence and extract that region specified by start and end positions for all sequences for bacteria and archaea seperately.
#   2. Merge all the truncated sequences together with the same header and the same order as seen in the original input FASTA file, and write it to file.
#        e.g.   ./InputData/ssu_all_r220_filtered.fna  ->  ./Output/ssu_all_r220_filtered_16s_v3.fna
# This is done for every region (e.g. V3, V3-V4, V6, V1-V3, etc.)
# ------------------------------

echo "Beginning extraction of 16S regions from alignments..."


# Define a helper function to compute the alignment range.
# This function scans the reference sequence (from the A2M file) and counts only aligned columns—
# i.e. those characters that are uppercase (A–Z) or the gap symbol (“-”).
# When the core count equals the requested start or end (from the .truncspec), it records the alignment index.
get_alignment_range() {
  local seq="$1"         # The full aligned reference sequence (a single string)
  local core_start="$2"  # The start position (from truncspec) in the reference's aligned columns
  local core_end="$3"    # The end position (from truncspec)
  local pos=0
  local core_count=0
  local aln_start=0
  local aln_end=0
  local char
  local seq_len=${#seq}
  for (( i=0; i<seq_len; i++ )); do
    char="${seq:$i:1}"
    # Count only uppercase letters and gap symbols ("-") as aligned columns.
    if [[ "$char" =~ [A-Z-] ]]; then
      (( core_count++ ))
      if [ "$core_count" -eq "$core_start" ]; then
        aln_start=$(( i+1 ))  # Convert to 1-indexed position
      fi
      if [ "$core_count" -eq "$core_end" ]; then
        aln_end=$(( i+1 ))
        break
      fi
    fi
  done
  echo "$aln_start $aln_end"
}

# Loop over each region (skip ARC_REF_SEQ_ID and BAC_REF_SEQ_ID)
for region in "${!region_specs[@]}"; do
  if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
    continue
  fi

  echo "Processing region: $region"
  # Parse parameters from the region spec string (e.g. "arc_start=512, arc_end=775, bac_start=346, bac_end=574, min_len=250, max_len=300")
  arc_start=$(echo "${region_specs[$region]}" | grep -oE 'arc_start=[0-9]+' | cut -d'=' -f2)
  arc_end=$(echo "${region_specs[$region]}" | grep -oE 'arc_end=[0-9]+' | cut -d'=' -f2)
  bac_start=$(echo "${region_specs[$region]}" | grep -oE 'bac_start=[0-9]+' | cut -d'=' -f2)
  bac_end=$(echo "${region_specs[$region]}" | grep -oE 'bac_end=[0-9]+' | cut -d'=' -f2)
  min_len=$(echo "${region_specs[$region]}" | grep -oE 'min_len=[0-9]+' | cut -d'=' -f2)
  max_len=$(echo "${region_specs[$region]}" | grep -oE 'max_len=[0-9]+' | cut -d'=' -f2)
  echo "  Bacteria region positions: $bac_start to $bac_end"
  echo "  Archaea region positions: $arc_start to $arc_end"
  echo "  Valid sequence length range: $min_len to $max_len bp"

  # ---- Bacterial reference and extraction ----
  # Find the bacterial reference record (its header should contain the $bac_ref substring)
  ref_bac_record=$(awk -v ref="$bac_ref" 'BEGIN{RS=">"; ORS=""} $0 ~ ref {print ">"$0; exit}' "${intermediates_dir}/align_bac_ssu.a2m")
  if [ -z "$ref_bac_record" ]; then
    echo "Error: Bacterial reference sequence not found for region $region."
    exit 1
  fi
  # Remove the header and newlines to get one continuous sequence string
  ref_bac_seq=$(echo "$ref_bac_record" | sed '1d' | tr -d '\n')
  # Compute the alignment index range for the bacterial reference using the helper function
  read bac_aln_start bac_aln_end < <(get_alignment_range "$ref_bac_seq" "$bac_start" "$bac_end")
  echo "  Bacterial alignment indices: $bac_aln_start to $bac_aln_end"

  # ---- Archaeal reference and extraction ----
  # Similarly, find the archaeal reference record (header contains $arc_ref)
  ref_arc_record=$(awk -v ref="$arc_ref" 'BEGIN{RS=">"; ORS=""} $0 ~ ref {print ">"$0; exit}' "${intermediates_dir}/align_arc_ssu.a2m")
  if [ -z "$ref_arc_record" ]; then
    echo "Error: Archaeal reference sequence not found for region $region."
    exit 1
  fi
  ref_arc_seq=$(echo "$ref_arc_record" | sed '1d' | tr -d '\n')
  read arc_aln_start arc_aln_end < <(get_alignment_range "$ref_arc_seq" "$arc_start" "$arc_end")
  echo "  Archaeal alignment indices: $arc_aln_start to $arc_aln_end"

  # ---- Extract truncated regions from bacterial and archaeal alignments ----
  # For each alignment file, we extract the substring from the computed alignment indices.
  # The extraction is done on the full record (FASTA format), then we post-process:
  #    - Remove any '-' and '.' characters
  #    - Convert the resulting sequence to uppercase

  echo "  Extracting truncated regions from bacterial and archaeal alignments..."

  # Extract bacterial truncated sequences into a temporary file
  truncated_bac_file="${intermediates_dir}/truncated_bac_${region}.fna"
  echo "  Extracting bacterial sequences to: $truncated_bac_file"
  awk -v start="$bac_aln_start" -v end="$bac_aln_end" 'BEGIN {RS=">"; ORS=""}
  NR > 1 {
    n = split($0, lines, "\n")
    header = lines[1]
    # Remove the index prefix (e.g., "0002|") from the header
    sub(/^[0-9]+\|/, "", header)
    seq = ""
    for(i=2; i<=n; i++){
      seq = seq lines[i]
    }
    sub_seq = substr(seq, start, end - start + 1)
    gsub(/[-.]/, "", sub_seq)
    # Convert to uppercase using the system tr command
    cmd = "echo " sub_seq " | tr \"[:lower:]\" \"[:upper:]\""
    cmd | getline sub_seq
    close(cmd)
    print ">" header "\n" sub_seq
  }' "${intermediates_dir}/align_bac_ssu.a2m" > "$truncated_bac_file"

  # Extract archaeal truncated sequences into a temporary file
  truncated_arc_file="${intermediates_dir}/truncated_arc_${region}.fna"
  echo "  Extracting archaeal sequences to: $truncated_arc_file"
  awk -v start="$arc_aln_start" -v end="$arc_aln_end" 'BEGIN {RS=">"; ORS=""}
  NR > 1 {
    n = split($0, lines, "\n")
    header = lines[1]
    # Remove the index prefix (e.g., "0002|") from the header
    sub(/^[0-9]+\|/, "", header)
    seq = ""
    for(i=2; i<=n; i++){
      seq = seq lines[i]
    }
    sub_seq = substr(seq, start, end - start + 1)
    gsub(/[-.]/, "", sub_seq)
    cmd = "echo " sub_seq " | tr \"[:lower:]\" \"[:upper:]\""
    cmd | getline sub_seq
    close(cmd)
    print ">" header "\n" sub_seq
  }' "${intermediates_dir}/align_arc_ssu.a2m" > "$truncated_arc_file"

  # ---- Merge truncated sequences in original FASTA order ----
  # The original input FASTA file has the desired header order.
  merged_output="${intermediates_dir}/02_reordered_truncated_${region}.fna"
  echo "  Merging bacterial and archaeal sequences in original FASTA order..."
  echo "  Writing merged output to: $merged_output"
  
    # Combine truncated bacterial and archaeal sequences into one temporary file
  joined_truncated_file="${intermediates_dir}/01_joined_truncated_${region}.fna"
  cat "$truncated_bac_file" "$truncated_arc_file" > "$joined_truncated_file"

  # Merge truncated sequences in original FASTA order.
  # We use the original input FASTA file to preserve the header order.
  # First, load all truncated sequences (from both bacteria and archaea)
  # into an associative array keyed by header (with the index prefix already removed).
  # Then, for each record in the input FASTA, we extract the header (without the ">" character)
  # and, if a corresponding truncated record exists, print the original header and its truncated sequence.
    # Merge truncated sequences in original FASTA order with debug output
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
        key = lines[1];  # Header (without the leading ">" because RS is ">" )
        seq = "";
        for (i = 2; i <= n; i++) {
          seq = seq lines[i];
        }
        truncated[key] = seq;  # Save the truncated sequence indexed by header
        count++;
      }
      next;
    }
    # Process the original FASTA file (second input file) to preserve header order
    {
      if ($0 != "") {
        n = split($0, lines, "\n");
        header = lines[1];  # This header should match one in truncated[]
        if (header in truncated) {
          print ">" header "\n" truncated[header] "\n";
          match_count++;
        }
      }
    }
  ' "$joined_truncated_file" "$input_fna" > "$merged_output"



  echo "  Region '$region' extraction complete. Output written to: $merged_output"
done

echo "Extraction completed for all regions."




# ------------------------------
# Filtering and Output Generation
# ------------------------------
# 1. For each region (e.g. V3, V3-V4), the truncated sequences are filtered 
#  - The truncated sequences found in ./Intermediates/reordered_truncated_{REGION}.fna are filtered
#  - Filters:
#    - Length: min_len-max_len+2*trunc_padding bp   (based on the .truncspec and trunc_padding)
#    - Ambiguous base content: 0%   (unless --no_filter_ambiguous is used)
#  - The truncated sequences that pass the filter are written to ./Intermediates/03_filtered_truncated_{REGION_NAME}.fna

# 2. The full sequences are filtered
#  - The full sequences from the original FASTA file (provided as $input_fna) are filtered
#  - Filters:
#    - Ambiguous base content: 0%   (unless --no_filter_ambiguous is used)
#  - The full sequences that pass the filter are written to ./Intermediates/filtered_full_seqs.fna
#  (regardless of whether --no_filter_ambiguous is used, ./Intermediates/filtered_full_seqs.fna is written)

# 3. Seqeuences are further filtered ensuring they passed for all regions
#  - If --no_require_all_regions flag is NOT given:
#    - For each sequence, we check if it is present in all regions (including the full sequences)
#      - This is done by looking in files: 
#        - ./Intermediates/filtered_full_seqs.fna
#        - ./Intermediates/03_filtered_truncated_{REGION_NAME}.fna
#    - If a sequence is present in all filtered files, it is written to the ./Output/*.fasta files
#    - This means that every ./Output/*.fasta files has the same set of sequences with the same headers in the same order as each other, but difference DNA sequences
#  - If --no_require_all_regions flag is given:
#    - These files are copied (and renamed) to ./Output/
#        - ./Intermediates/filtered_full_seqs.fna -> ./Output/FULL_seqs.fasta
#        - ./Intermediates/03_filtered_truncated_{REGION_NAME}.fna -> ./Output/{REGION}_seqs.fasta
#  - All .fasta output files are written to ./Output/
#    - ./Output/{REGION}_seqs.fasta
#    - ./Output/FULL_seqs.fasta
#  (regardless of whether --no_require_all_regions is used, all output files are written)

# 4. The about.txt file is written
#  - The ./Output/about.txt file is written which summarizes the entire processing run of the 16S rRNA region extraction.
#  - The file contains the following information:
#    - Lists the key input parameters used by the script:
#      • Input FASTA file: The original sequence file.
#      • Bacteria HMM and Archaea HMM: The Hidden Markov Model files used for sequence alignment.
#      • Truncation specification file: Defines the regions to extract with their start/end positions
#        and length boundaries.
#      • Option flags: Indicates if alignment was skipped, ambiguous filtering was disabled,
#        whether sequences from all regions are required, if intermediate files were removed,
#        and the value of truncation padding.
#    - Reports the reference sequence IDs for Archaea and Bacteria.
#    - Iterates through each defined region (e.g., V3, V3-V4, etc.) and prints its extraction
#      parameters such as start and end positions, and the expected minimum and maximum lengths.
#    - Displays the total number of sequences from the original input FASTA.
#    - Shows how many sequences passed all filters overall, along with the corresponding percentage.
#    - For each region, prints:
#      • The number of sequences that passed the filters.
#      • The number (and percentage) of sequences that were filtered out.
#    - Lists the specific filtering criteria for each region:
#      • The length criteria (min-max bp).
#      • Whether the ambiguous base content filter was applied (<1% ambiguous bases).
#      • Whether sequences were required to be present in all regions.
#    - Logs the date and time when the processing was completed.
# ------------------------------------------------------------------------------

echo "Starting filtering and output generation process..."

# Define an output directory for final output files
output_dir="./Output"
mkdir -p "$output_dir"
echo "Created output directory: $output_dir"

#############################################
# 1. Filter truncated sequences per region
#############################################
for region in "${!region_specs[@]}"; do
  # Skip the ARC_REF_SEQ_ID and BAC_REF_SEQ_ID entries
  if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
    continue
  fi
  
  echo "Filtering truncated sequences for region: $region"
  # Extract min and max lengths from the region specification
  min_len=$(echo "${region_specs[$region]}" | grep -oE 'min_len=[0-9]+' | cut -d'=' -f2)
  max_len=$(echo "${region_specs[$region]}" | grep -oE 'max_len=[0-9]+' | cut -d'=' -f2)
  # Compute effective length bounds (each side gets trunc_padding)
  lower_bound=$(( min_len + 2 * trunc_padding ))
  upper_bound=$(( max_len + 2 * trunc_padding ))
  
  input_file="${intermediates_dir}/02_reordered_truncated_${region}.fna"
  output_file="${intermediates_dir}/03_filtered_truncated_${region}.fna"
  
  if [ "$no_filter_ambiguous" = true ]; then
    # Only apply length filtering when ambiguous base filtering is disabled
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
  echo "  Filtered truncated sequences written to: $output_file"
done

#############################################
# 2. Filter full sequences from the input FASTA
#############################################
echo "Filtering full sequences from input FASTA..."
full_filtered="${intermediates_dir}/filtered_full_seqs.fna"
if [ "$no_filter_ambiguous" = true ]; then
  # No ambiguous filtering: simply copy the input FASTA
  cp "$input_fna" "$full_filtered"
else
  awk 'BEGIN { RS=">"; ORS="" }
  {
    if($0 != ""){
      n = split($0, lines, "\n")
      header = lines[1]
      seq = ""
      for(i=2; i<=n; i++){
        seq = seq lines[i]
      }
      # Only print records with sequences of A, T, C, G only
      if(seq ~ /^[ATCG]+$/){
        print ">" header "\n" seq "\n"
      }
    }
  }' "$input_fna" > "$full_filtered"
fi
echo "  Filtered full sequences written to: $full_filtered"

#############################################
# 3. Cross-region filtering (or direct copy if --no_require_all_regions)
#############################################
if [ "$no_require_all_regions" = true ]; then
  echo "Skipping cross-region filtering as --no_require_all_regions flag is set."
  # Copy the filtered full sequences and each region’s truncated sequences to the Output directory
  cp "$full_filtered" "$output_dir/FULL_seqs.fasta"
  for region in "${!region_specs[@]}"; do
    if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
      continue
    fi
    region_filtered="${intermediates_dir}/03_filtered_truncated_${region}.fna"
    cp "$region_filtered" "$output_dir/${region}_seqs.fasta"
    echo "  Copied filtered truncated sequences for region $region to: $output_dir/${region}_seqs.fasta"
  done
else
  echo "Performing cross-region filtering to retain only sequences present in all filtered files..."
  # Start by extracting headers from the filtered full sequences file
  temp_common="${intermediates_dir}/common_headers.txt"
  grep '^>' "$full_filtered" | sed 's/^>//' > "$temp_common"
  
  # For each region, intersect its headers with the common header list
  for region in "${!region_specs[@]}"; do
    if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
      continue
    fi
    region_file="${intermediates_dir}/03_filtered_truncated_${region}.fna"
    temp_headers="${intermediates_dir}/headers_${region}.txt"
    grep '^>' "$region_file" | sed 's/^>//' > "$temp_headers"
    intersect="${intermediates_dir}/common_tmp.txt"
    grep -Fxf "$temp_common" "$temp_headers" > "$intersect"
    mv "$intersect" "$temp_common"
    rm "$temp_headers"
  done
  echo "  Common headers count: $(wc -l < "$temp_common")"
  
  # Define a helper function to filter a FASTA file by common headers
  filter_by_common_headers() {
    local infile="$1"
    local outfile="$2"
    awk -v common="$temp_common" 'BEGIN{
      # Read common headers into an array
      while((getline line < common) > 0){
        headers[line]=1
      }
      RS=">"; ORS=""
    }
    {
      if($0 != ""){
        n = split($0, lines, "\n")
        header = lines[1]
        if(header in headers){
          printf(">%s\n", $0)
        }
      }
    }' "$infile" > "$outfile"
  }
  
  # Filter the full sequences and each region’s truncated sequences using the common headers
  filter_by_common_headers "$full_filtered" "$output_dir/FULL_seqs.fasta"
  echo "  Full sequences after cross-region filtering written to: $output_dir/FULL_seqs.fasta"
  for region in "${!region_specs[@]}"; do
    if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
      continue
    fi
    region_file="${intermediates_dir}/03_filtered_truncated_${region}.fna"
    output_region="$output_dir/${region}_seqs.fasta"
    filter_by_common_headers "$region_file" "$output_region"
    echo "  Filtered truncated sequences for region $region written to: $output_region"
  done
  
  # Clean up temporary common headers file
  rm "$temp_common"
fi

#############################################
# 4. Generate about.txt summary file
#############################################
echo "Generating about.txt summary file..."
{
  echo "16S rRNA Region Extraction Summary"
  echo "=================================="
  echo ""
  echo "Input Parameters:"
  echo "  Input FASTA file: $input_fna"
  echo "  Bacteria HMM: $bac_hmm"
  echo "  Archaea HMM: $arc_hmm"
  echo "  Truncation specification file: $trunc_spec_file"
  echo "  Options:"
  echo "    Skip alignment: $skip_align"
  if [ "$no_filter_ambiguous" = true ]; then
    echo "    Filter ambiguous bases: Disabled"
  else
    echo "    Filter ambiguous bases: Enabled"
  fi
  if [ "$no_require_all_regions" = true ]; then
    echo "    Require all regions: Not required"
  else
    echo "    Require all regions: Required"
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
    echo "    Expected sequence length: $effective_min to $effective_max bp"
    if [ "$no_filter_ambiguous" = true ]; then
      echo "    Ambiguous base filter: Not applied"
    else
      echo "    Ambiguous base filter: Applied (0% ambiguous bases allowed)"
    fi
  done
  echo ""
  total_seqs=$(grep -c "^>" "$input_fna")
  echo "Total sequences in input FASTA: $total_seqs"
  if [ "$no_require_all_regions" = true ]; then
    passed_seqs=$(grep -c "^>" "$full_filtered")
  else
    passed_seqs=$(grep -c "^>" "$output_dir/FULL_seqs.fasta")
  fi
  percentage=$(awk "BEGIN {printf \"%.2f\", ($passed_seqs/$total_seqs*100)}")
  echo "Sequences passed all filters: $passed_seqs ($percentage%)"
  echo ""
  for region in "${!region_specs[@]}"; do
    if [[ "$region" == "ARC_REF_SEQ_ID" || "$region" == "BAC_REF_SEQ_ID" ]]; then
      continue
    fi
    if [ "$no_require_all_regions" = true ]; then
      region_file="${intermediates_dir}/03_filtered_truncated_${region}.fna"
    else
      region_file="$output_dir/${region}_seqs.fasta"
    fi
    region_passed=$(grep -c "^>" "$region_file")
    filtered_out=$(( total_seqs - region_passed ))
    perc_filtered=$(awk "BEGIN {printf \"%.2f\", ($filtered_out/$total_seqs*100)}")
    echo "Region $region: $region_passed passed, $filtered_out filtered out ($perc_filtered% filtered)"
  done
  echo ""
  echo "Processing completed on: $(date)"
} > "$output_dir/about.txt"
echo "about.txt written to: $output_dir/about.txt"






# Clean up intermediate files if requested
if [ "$rm_intermediates" = true ]; then
  echo "Cleaning up intermediate files..."
  rm -rf "${intermediates_dir}"
  echo "Intermediate files removed."
else
  echo "Intermediate files kept in: ${intermediates_dir}"
  echo "Use --rm_intermediates flag to remove these files automatically."
fi

# Calculate and display execution time
end_time=$(date +%s)
execution_time=$((end_time - start_time))
minutes=$((execution_time / 60))
seconds=$((execution_time % 60))

echo "Done. All output files have been written to: $output_dir"
echo "Total execution time: ${minutes} minutes and ${seconds} seconds."
