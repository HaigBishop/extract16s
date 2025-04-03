
# extract16s

A tool for extracting variable regions from 16S rRNA sequences.

## Overview

`extract16s` is a command-line tool that uses Hidden Markov Models (HMMs) to extract specific variable regions from 16S rRNA gene sequences. The tool aligns input sequences to reference HMMs and then extracts regions of interest defined in a truncation specification file (.truncspec).

Key features:
- Handles bacterial and archaeal 16S rRNA sequences seperately
- Enabled use of reference sequences or truncation positioning
- Uses HMMs for alignment and region extraction
- Filters sequences by length and ambiguous base content
- Extracts multiple variable regions (V3, V4, V3-V4, etc.) in a single run

## Requirements

- **HMMER suite**: Provides `hmmalign` and `esl-reformat`. [User Guide](http://eddylab.org/software/hmmer/Userguide.pdf)
- **Unix tools**: awk, bc, grep, sed, tr

## Installation

```bash
git clone https://github.com/HaigBishop/extract16s.git
cd extract16s
```

## Directory Structure

The tool expects a specific directory structure:

```
extract16s/
├── InputData/
│   ├── arc_16s.hmm                # Archaeal 16S HMM file
│   ├── bac_16s.hmm                # Bacterial 16S HMM file
│   ├── ssu_all_filtered.fna       # Input FASTA file
│   └── trunc_spec_v3_v3-v4.truncspec  # Truncation specification file
│
├── Scripts/
│   └── extract16s.sh              # Main script
│
├── Intermediates/                 # Created during execution
│   └── ...                        # Various intermediate files
│
└── Output/                        # Final output directory
    ├── FULL_seqs.fasta           # Filtered full-length sequences
    ├── V3_seqs.fasta             # Extracted V3 region sequences
    ├── V3-V4_seqs.fasta          # Extracted V3-V4 region sequences
    └── about_extraction.txt                 # Summary of processing details
```

## Usage

```bash
bash ./Scripts/extract16s.sh <input_fna> <bac_hmm> <arc_hmm> <trunc_spec_file> [options]
```

### Required Arguments

- `<input_fna>`: Input FASTA file containing 16S rRNA sequences
- `<bac_hmm>`: Bacterial 16S rRNA HMM file
- `<arc_hmm>`: Archaeal 16S rRNA HMM file
- `<trunc_spec_file>`: Truncation specification file (.truncspec)

### Options

- `--skip_align`: Skip alignment process (use existing alignments in Intermediates/)
- `--no_filter_ambiguous`: Skip filtering by ambiguous base content
- `--no_require_all_regions`: Don't require sequences to be successfully extracted for all regions
- `--rm_intermediates`: Remove intermediate files after processing
- `--trunc_padding N`: Add N bases of padding to each side of extracted regions (default: 0)

### Example

```bash
bash ./Scripts/extract16s.sh \
     ./InputData/ssu_all_filtered.fna \
     ./InputData/bac_16s.hmm \
     ./InputData/arc_16s.hmm \
     ./InputData/trunc_spec_v3_v3-v4.truncspec \
     --trunc_padding 15 --rm_intermediates
```

## Input Files

### FASTA File Format

The input FASTA file should contain 16S rRNA gene sequences with headers that indicate whether they are bacterial or archaeal (containing either `d__Bacteria` or `d__Archaea`).

Example:
```
>GB_GCA_000488855.1~KI535272.1 d__Bacteria;p__Bacillota_A;c__Clostridia;o__Peptostreptococcales;f__Anaerovoracaceae;g__Gallibacter;s__Gallibacter brachus [location=411..1920] [ssu_len=1510] [contig_len=1959]
GGCAGCAGTGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAA...
>GB_GCA_000494145.1~AWOE01000013.1 d__Archaea;p__Thermoproteota;c__Nitrososphaeria_A;o__Caldarchaeales;f__Calditenuaceae;g__Calditenuis;s__Calditenuis sp000494145 [location=10546..12034] [ssu_len=1489] [contig_len=32668]
CGCCCGTAGCCGGCCCGGTGTGTCCCTCGTTAAATCCACGGGCTTAACCCGTGGGCTGC...
```

### Truncation Specification File (.truncspec)

The truncation specification file defines the reference sequence IDs and the regions to be extracted:

```
# 16S Region Truncation Specification File (.truncspec)
ARC_REF_SEQ_ID: RS_GCF_000025685.1~NC_013967.1
BAC_REF_SEQ_ID: RS_GCF_000194175.1~NZ_AEZI02000002.1
V4: arc_start=471, arc_end=724, bac_start=530, bac_end=782, min_len=245, max_len=260
V3-V4: arc_start=336, arc_end=722, bac_start=354, bac_end=780, min_len=385, max_len=445
```

This file must include:
- `ARC_REF_SEQ_ID`: Reference sequence ID for archaea
- `BAC_REF_SEQ_ID`: Reference sequence ID for bacteria
- One or more region definitions, each with:
  - Region name (e.g., V4, V3-V4)
  - Start and end positions for archaeal reference (arc_start, arc_end)
  - Start and end positions for bacterial reference (bac_start, bac_end)
  - Minimum and maximum expected sequence lengths (min_len, max_len)

### HMM Files

The tool requires two HMM files for alignment:

- `InputData/bac_16s.hmm`: Bacterial 16S rRNA HMM
- `InputData/arc_16s.hmm`: Archaeal 16S rRNA HMM

These HMM files are taken from the [barrnap](https://github.com/tseemann/barrnap/) tool's [database](https://github.com/tseemann/barrnap/blob/master/db/). They contain profile Hidden Markov Models specifically trained for bacterial and archaeal 16S rRNA genes. The files have been modified to remove non-16S models (5S, 23S, etc.).

## Output Files

The tool generates several output files in the `Output/` directory:

- `FULL_seqs.fasta`: Full-length sequences that passed filters
- `{REGION}_seqs.fasta`: Extracted sequences for each region (e.g., V3_seqs.fasta, V3-V4_seqs.fasta)
- `about_extraction.txt`: Summary of processing details, including:
  - Input parameters
  - Reference sequence IDs
  - Region extraction parameters
  - Filtering statistics
  - Processing timestamps

## Processing Steps

1. Sequences are split into bacterial and archaeal sets
2. Each set is aligned to its corresponding HMM
3. Alignments are converted to A2M format
4. Variable regions are extracted based on the truncation specification
5. Sequences are filtered by:
   - Length (based on min_len, max_len, and padding)
   - Ambiguous base content (unless --no_filter_ambiguous is used)
6. Optionally, sequences absent in any region are filtered out
7. Final files are written to the Output directory


## License

This tool is provided under the MIT License.

## Contact

For questions, issues, or suggestions, please contact:
- Email: haigvbishop@gmail.com

---
Created by Haig Bishop, 2025
