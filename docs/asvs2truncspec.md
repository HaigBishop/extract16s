# The asvs2truncspec Tool

The `asvs2truncspec` tool creates `.truncspec` files from ASV (Amplicon Sequence Variant) sequences. A `.truncspec` file specifies exact regions to truncate full length 16S genes when using `extract16s.sh`.

This tool takes a directory of ASV sequences from any number of different 16S regions and creates a `.truncspec` file that specifies regions matching the provided ASV sequences. The tool is implemented as a combination of Bash and Python scripts: `asvs2truncspec_01_prep.py`, `asvs2truncspec_02_hmmer.sh`, and `asvs2truncspec_03_agg.py`.


---


## Input Data

### Full length 16S GTDB database
The tool requires the full length 16S GTDB database, for example `ssu_all_r226.fna`. The path to this file is set as a configuration constant `DB_PATH` in the Python scripts. The format is standard FASTA with GTDB headers:
```fasta
>RS_GCF_002194975.1~NZ_NHWV01000143.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli [location=150..1687] [ssu_len=1538] [contig_len=1698]
TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGAAGCTTGCTTCTTTGCTGACGAGTGGC
>RS_GCF_038435785.1~NZ_JARMYW010000081.1 d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli [location=5..1448] [ssu_len=1444] [contig_len=1463]
ACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAG
```

### ASV sequences and dataset metadata
The tool requires ASV sequences from any number of different 16S regions, organized in a directory specified by `ASV_DIR`. Each FASTA file contains ASV sequences from one or more datasets:
```bash
target_asvs/
├── AGP.fna
├── CMR-Trunc_DA.fna
├── CMR-Trunc_SA.fna
├── CMR-Trunc_SP.fna
├── CMR-Trunc_ST.fna
├── CMR-Untrunc_DA.fna
├── CMR-Untrunc_SA.fna
├── CMR-Untrunc_SP.fna
├── CMR-Untrunc_ST.fna
├── HMC.fna
└── SAR.fna
```

Example ASV sequences:
```fasta
>SAR_ASV_1
TACGGAAGGTCCGGGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGCCGGAGATTAAGCGTGTTGTGAAATGTAGACGCTCAACGTCTGCACTGCAGCGCGAACTGGTTTCCTTGAGTACGCACAAAGTGGGCGGAATTCGTGG
>SAR_ASV_2
TACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTCTGTCAAGTCGGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATTCGAAACTGGCAGGCTAGAGTCTTGTAGAGGGGGGTAGAATTCCAGG
```

Optionally, TSV metadata files can specify which dataset each ASV belongs to. This is important because the tool requires that ASVs from the same dataset are all from the exact same 16S region. A FASTA file's corresponding TSV file must have the same name as the FASTA file with `_datasets.tsv` appended instead of `.fna`:
```bash
target_asvs/
├── HMC_datasets.tsv
├── HMC.fna
├── CMR-Untrunc_SP_datasets.tsv
├── CMR-Untrunc_SP.fna
└── ...
```

The TSV files must include the columns `ASV_ID` and `Dataset_ID` with rows for every ASV in the corresponding FASTA file. An optional `Region` column can specify the 16S region (e.g. "V4", "V3-V4") for each ASV:
```tsv
ASV_ID	Dataset_ID	Region
CMR_SP_ASV_1	16S_68_Girdhar	V4
CMR_SP_ASV_2	16S_68_Girdhar	V4
CMR_SP_ASV_3	16S_68_Girdhar	V4
CMR_SP_ASV_4	16S_1211_Milletich	V3-V4
CMR_SP_ASV_5	16S_1211_Milletich	V3-V4
```

If no TSV file is provided for a given FASTA file, the tool assumes that all ASVs in the FASTA file are from the same dataset. If the `Region` column is not present (or is not provided for some ASVs), the region is set to "UNK" (unknown) in the dataset manifest.

### 16S HMMs (Hidden Markov Models)
The tool requires two HMM files for alignment:
- `InputData/bac_16s.hmm`: Bacterial 16S rRNA HMM
- `InputData/arc_16s.hmm`: Archaeal 16S rRNA HMM

These HMM files are taken from the [barrnap](https://github.com/tseemann/barrnap/) tool's [database](https://github.com/tseemann/barrnap/blob/master/db/). They contain profile Hidden Markov Models specifically trained for bacterial and archaeal 16S rRNA genes. The files have been modified to remove non-16S models (5S, 23S, etc.).

### Reference Sequences
The tool requires archaeal and bacterial reference sequences from the GTDB database. All final position coordinates are relative to these references (rather than using alignment coordinates). The reference sequence IDs must be provided as configuration constants (`ARC_REF_SEQ_ID` and `BAC_REF_SEQ_ID`), and they must be found in the DB_PATH file. IDs take the form `RS_GCF_022846175.1~NZ_AP025587.1-#2`. The tool extracts the first sequence whose header contains the ID string.


---


## Configuration

These configuration values control the tool's behavior and are not specific to any one step in the pipeline. They are set as global constants at the top of the Python scripts.

To configure the tool, provide paths to the input data and output data:
```python
DB_PATH = "/mnt/secondary/micro16s_dbs/16s_databases/ssu_all_r226.fna"
ASV_DIR = "/home/haig/Repos/micro16s/extract16s/target_asvs"
TRUNCSPEC_OUT_PATH = "/home/haig/Repos/micro16s/extract16s/asvs2truncspec_out/new.truncspec"
INFO_OUT_DIR = "/home/haig/Repos/micro16s/extract16s/asvs2truncspec_out/"

# HMMs (used by the Step 2 bash script)
BAC_HMM_PATH = "InputData/bac_16s.hmm"
ARC_HMM_PATH = "InputData/arc_16s.hmm"

# Reference sequences
ARC_REF_SEQ_ID = "RS_GCF_022846175.1~NZ_AP025587.1-#2"
BAC_REF_SEQ_ID = "RS_GCF_030545895.1~NZ_JAUOMX010000042.1"

# Region length buffers
MIN_LEN_BUFFER = 50
MAX_LEN_BUFFER = 50

# Optional cross-domain bootstrapping (Step 3)
USE_CROSS_DOMAIN_BOOTSTRAPPING = True
CROSS_DOMAIN_BOOTSTRAPPING_MAX_DISTANCE = 90

# Optional region redundancy minimisation (Step 3)
USE_REGION_REDUNDANCY_MINIMISATION = True
REGION_REDUNDANCY_MINIMISATION_MAX_DISTANCE = 20
REGION_REDUNDANCY_MINIMISATION_NAME_JOINER = "+"

# Optional region renaming (Step 3)
USE_REGION_RENAMING = True
REGION_RENAMING_METHOD = "uniform"  # "uniform", "dataset", or "region"

# Intermediate directory structure (all derived from INFO_OUT_DIR)
INTER_DIR = f"{INFO_OUT_DIR}/intermediates"
STAGED_DIR = f"{INTER_DIR}/01_staged_inputs"
HMMER_OUT_DIR = f"{INTER_DIR}/02_hmmer"
RESULTS_DIR = f"{INTER_DIR}/03_results"

# Dataset metadata file conventions
DATASET_TSV_SUFFIX = "_datasets.tsv"
DATASET_TSV_ASV_COL = "ASV_ID"
DATASET_TSV_DATASET_COL = "Dataset_ID"
DATASET_TSV_REGION_COL = "Region"  # optional column
UNKNOWN_REGION = "UNK"  # value used when Region is not available

# Step 2 output size guard
NHMMER_REPORT_MAX_EVALUE = 1e-3

# Robust hit filtering / consensus calling
ARC_HMMER_MAX_EVALUE = 1e-4
BAC_HMMER_MAX_EVALUE = 1e-5
MIN_HIT_COVERAGE_FRAC = 0.6
ARC_MIN_HITS_PER_DATASET = 5
BAC_MIN_HITS_PER_DATASET = 5
OUTLIER_TOL_BP = 30

# Number of threads to use for HMMER searches
THREADS = 18
```


---


## Output Data

### `.truncspec` file
The tool outputs a `.truncspec` file that specifies the regions to truncate for each dataset. Lines starting with `#` are comments and are ignored by `extract16s.sh`. The format is:
```truncspec
# 16S Region Truncation Specification File (.truncspec)
# Generated by asvs2truncspec using 267 datasets to derive 31 regions
# This file specifies regions to be extracted from 16S alignments using `extract16s.sh`

ARC_REF_SEQ_ID: RS_GCF_022846175.1~NZ_AP025587.1-#2
BAC_REF_SEQ_ID: RS_GCF_030545895.1~NZ_JAUOMX010000042.1
AGP: arc_start=515, arc_end=716, bac_start=412, bac_end=637, min_len=152, max_len=276
SAR: arc_start=390, arc_end=615, bac_start=569, bac_end=794, min_len=176, max_len=276
HMC_Dataset_001: arc_start=515, arc_end=716, bac_start=412, bac_end=637, min_len=152, max_len=276
HMC_Dataset_002: arc_start=390, arc_end=615, bac_start=569, bac_end=794, min_len=176, max_len=276
CMR-Untrunc_SP_Girdhar: arc_start=515, arc_end=716, bac_start=412, bac_end=637, min_len=152, max_len=276
CMR-Untrunc_SP_Milletich: arc_start=390, arc_end=615, bac_start=569, bac_end=794, min_len=176, max_len=276
...
```

Notice how the regions correspond to datasets (not input files). For instance, there is exactly one AGP region because the `AGP.fna` file had no supplied dataset metadata file, so all ASVs in the AGP file are from the same dataset. For CMR-Untrunc_SP, there are two regions because there are two datasets in the `CMR-Untrunc_SP.fna` file, as specified in the `CMR-Untrunc_SP_datasets.tsv` file.

If cross-domain bootstrapping is enabled, two `.truncspec` files are written:
- The bootstrapped file at `TRUNCSPEC_OUT_PATH`
- An unbootstrapped copy with `_no-bootstrapping` inserted before the extension (e.g. `new_no-bootstrapping.truncspec`)

If region redundancy minimisation is enabled, an additional `.truncspec` file is written with `_no-redundancy-min` inserted before the extension. This is the pre-merge version used for comparison.

If region renaming is enabled, an additional `.truncspec` file is written with `_no-renaming` inserted before the extension. This is the pre-rename version used for comparison. Additionally, a `region_renamings.tsv` file is written to `INFO_OUT_DIR/intermediates/03_results/` which maps old names to new names.

Combined output behavior:
- All three features enabled: 4 files (`_no-bootstrapping`, `_no-redundancy-min`, `_no-renaming`, final).
- Redundancy minimisation + region renaming enabled: 3 files (`_no-redundancy-min`, `_no-renaming`, final).
- Redundancy minimisation + bootstrapping enabled: 3 files (`_no-bootstrapping`, `_no-redundancy-min`, final).
- Region renaming + bootstrapping enabled: 3 files (`_no-bootstrapping`, `_no-renaming`, final).
- Only one feature enabled: 2 files (pre-feature file, final).
- No optional features enabled: 1 file (final only).

### Final region plots
If at least one region has complete coordinates for both domains, Step 3 writes two plot images directly to `INFO_OUT_DIR`:
- `final_regions_arc.png`
- `final_regions_bac.png`

These plots show the final truncspec regions (only those with both domains intact) across the full archaeal or bacterial reference sequence. The x-axis is reference base position with ticks every 50 bp. The y-axis is the region name list (long names truncated to 12 characters plus `...`), and each region is drawn as a simple blue rectangle.

#### Region Starts and Ends
The primary product of the `asvs2truncspec` tool is start and end coordinates in the `.truncspec` file. Separately for Archaeal and Bacterial reference sequences, the tool finds the start and end coordinates of the region that matches the ASV sequences best. The `extract16s.sh` script uses these coordinates to extract the region from the 16S alignment by converting the sequence coordinates to alignment coordinates and effectively truncating the whole alignment.

#### Maximum and minimum region lengths
The `extract16s.sh` script requires minimum and maximum region lengths to be specified. These act as a filter to remove extraction of regions that are unreasonably short or long. The `asvs2truncspec` tool calculates `min_len` and `max_len` based on the following logic:
- If both archaeal and bacterial regions are found:
    - `min_len = min(arc_end - arc_start + 1, bac_end - bac_start + 1) - MIN_LEN_BUFFER`
    - `max_len = max(arc_end - arc_start + 1, bac_end - bac_start + 1) + MAX_LEN_BUFFER`
- If only one domain's region is found:
    - `min_len = (domain_end - domain_start + 1) - MIN_LEN_BUFFER`
    - `max_len = (domain_end - domain_start + 1) + MAX_LEN_BUFFER`

The `MIN_LEN_BUFFER` and `MAX_LEN_BUFFER` configuration constants can be set to adjust the forgiveness for variation in region lengths.

#### Failure to find region coordinates
When the tool cannot confidently call region coordinates for a given dataset/domain, Step 3 records `NA` in the QC tables (e.g. `dataset_region_calls.tsv` / `dataset_region_qc.tsv`). This most commonly happens when:
- the dataset contains few/no ASVs from that domain (e.g. bacteria-dominated datasets often have too few archaeal hits), or
- HMMER hits are too weak/inconsistent to form a consensus region (after filtering/outlier logic).

Because `extract16s.sh` expects numeric coordinates, `NA` values cannot be used as input to `extract16s.sh`.

Therefore, when writing the final `.truncspec` file, any region with incomplete coordinates is still written for visibility, but it is **commented out** (prefixed with `#`) so `extract16s.sh` ignores it.

Example (commented-out incomplete regions):
```truncspec
ARC_REF_SEQ_ID: RS_GCF_022846175.1~NZ_AP025587.1-#2
BAC_REF_SEQ_ID: RS_GCF_030545895.1~NZ_JAUOMX010000042.1
#HMC_Dataset_001: arc_start=NA, arc_end=NA, bac_start=412, bac_end=637, min_len=152, max_len=276
#HMC_Dataset_002: arc_start=NA, arc_end=NA, bac_start=NA, bac_end=NA, min_len=NA, max_len=NA
```


### Logging/info/plots
The tool also outputs logging/info to the directory specified by `INFO_OUT_DIR`.


---


## Pipeline

This pipeline keeps HMMER calls in bash (since they are CLI tools), and keeps everything else in Python. The overall process is:
1) Stage ASVs and reference sequences into a predictable intermediate directory, with explicit dataset IDs.
2) Run two HMMER searches (archaea + bacteria) to get **model coordinates** for each ASV, plus two HMM-to-reference alignments for the reference sequences.
3) Convert model coordinates to **reference-sequence coordinates**, then compute a robust consensus start/end for each dataset and write the `.truncspec`.

### Pipeline Scripts

The `asvs2truncspec` tool is implemented as three numbered scripts so that each step is understandable and rerunnable:
- `Scripts/asvs2truncspec_01_prep.py` (Python): staging/validation + write intermediate inputs
- `Scripts/asvs2truncspec_02_hmmer.sh` (bash): run HMMER tools and write raw outputs
- `Scripts/asvs2truncspec_03_agg.py` (Python): parse HMMER outputs + compute consensus regions + write `.truncspec`

Design choices:
- Bash scripts expose CLI arguments (so they are easy to call/reuse).
- Python scripts keep all configuration options as global constants at the top of the file.

The pipeline creates a dedicated working directory under `INFO_OUT_DIR`:
```bash
asvs2truncspec_out/
├── intermediates/
│   ├── 01_staged_inputs/
│   ├── 02_hmmer/
│   └── 03_results/
├── logs/
├── final_regions_arc.png
└── final_regions_bac.png
```

The general rule is: Step 1 creates `01_staged_inputs/`, Step 2 reads that directory and creates `02_hmmer/`, and Step 3 reads both and creates `03_results/` plus the final `.truncspec`.

### Coordinate Conventions
To avoid off-by-one confusion, this pipeline uses the following conventions:
- **`.truncspec` coordinates (`arc_start`/`arc_end`/`bac_start`/`bac_end`)**: **1-based, inclusive**, referring to nucleotide positions in the *unaligned* reference sequences (this matches what `extract16s.sh` expects).
- **HMMER model coordinates (`hmmfrom`/`hmmto`)**: treated as **1-based, inclusive** as reported in `nhmmer --tblout` output.
- **Mapping tables** (e.g. `arc_model_to_ref.tsv` / `bac_model_to_ref.tsv`): store **model_pos (1-based)** → **ref_pos (1-based)** so they can be used directly with tblout positions and `.truncspec` output.
- **Internal Python indexing**: free to use 0-based indices for strings/arrays as normal, but any persisted coordinates (TSVs and `.truncspec`) are written as **1-based, inclusive** for simplicity and consistency.
- **Region lengths**: always use inclusive coordinates when computing lengths: `len = end - start + 1`.

### Processing Step 1 (asvs2truncspec_01_prep.py)

This step prepares the inputs into a clean, predictable structure, so later steps do not have to keep re-discovering metadata.

**Inputs:**
- `DB_PATH` (full length 16S GTDB FASTA, e.g. `ssu_all_r226.fna`)
- `ARC_REF_SEQ_ID`, `BAC_REF_SEQ_ID` (must be present in `DB_PATH`)
- `ASV_DIR` containing `*.fna` and optional `*_datasets.tsv`

**Validation rules:**
- Missing ASVs in a TSV file (present in FASTA but absent from TSV) is a hard fail.
- ASVs present in TSV but absent from FASTA are allowed (no warning).
- Duplicate ASV IDs are a hard fail.
- Dataset/ASV/source IDs must not contain `|`, `=`, whitespace, or non-printing characters; violations are a hard fail.

**Outputs (all under `INFO_OUT_DIR/intermediates/01_staged_inputs/`):**

1) **Reference sequences extracted from the DB**
```bash
ref_arc.fna
ref_bac.fna
```

2) **One combined FASTA of all ASVs**

This file is used for the HMMER searches so we only run each HMM once:
```bash
all_asvs.fna
```
Headers are rewritten to be unambiguous and machine-parseable:
```fasta
>dataset=HMC_Dataset_001|asv=ASV105105|src=HMC.fna
ACGT...
```

3) **Dataset manifest**

This is the source of truth for which ASV belongs to which dataset, including the 16S region metadata:
```bash
dataset_manifest.tsv
```
```tsv
dataset_id	asv_id	source_fna	region
HMC_Dataset_001	ASV105105	HMC.fna	V4
HMC_Dataset_001	ASV105106	HMC.fna	V4
CMR-Untrunc_SP_Girdhar	CMR_SP_ASV_1	CMR-Untrunc_SP.fna	V3-V4
CMR-Untrunc_SP_Milletich	CMR_SP_ASV_4	CMR-Untrunc_SP.fna	UNK
```

If there is no `*_datasets.tsv` for a FASTA file, Step 1 assumes that all ASVs in that FASTA file belong to a single dataset, and uses the FASTA base name as the dataset name (e.g. `AGP` for `AGP.fna`). If a `*_datasets.tsv` file is present, each dataset is named `{fna_base}_{Dataset_ID}`.

The `region` column contains the 16S region (e.g. "V4", "V3-V4") from the input TSV's `Region` column if available, otherwise "UNK" (unknown).

4) **Staging information**

A summary of the staged data, reference sequences, and the configuration used:
```bash
about_staging.txt
```

**Running Step 1:**
```bash
python Scripts/asvs2truncspec_01_prep.py
```

### Processing Step 2 (asvs2truncspec_02_hmmer.sh)

This step runs the external HMMER tools and writes their raw outputs. The key output from the ASV searches is **HMM model coordinates** (`hmmfrom`/`hmmto`) for each ASV hit, because model coordinates are stable across taxa and easy to convert into reference-sequence coordinates later.

**Inputs:**
- `InputData/arc_16s.hmm`
- `InputData/bac_16s.hmm`
- `INFO_OUT_DIR/intermediates/01_staged_inputs/ref_arc.fna`
- `INFO_OUT_DIR/intermediates/01_staged_inputs/ref_bac.fna`
- `INFO_OUT_DIR/intermediates/01_staged_inputs/all_asvs.fna`

**Outputs (all under `INFO_OUT_DIR/intermediates/02_hmmer/`):**

1) **Reference-to-HMM alignments**

Alignments of the two reference sequences to their HMMs so that Step 3 can build a mapping from HMM model positions → reference sequence positions:
```bash
arc_ref.sto
bac_ref.sto
```
(These are Stockholm outputs from `hmmalign`.)

2) **ASV-to-HMM search results**

Searches of the ASVs against each HMM with `nhmmer`. Uses `--tblout` output because it includes stable, parseable alignment/model coordinates (`hmmfrom`/`hmmto`, `alifrom`/`alito`, `envfrom`/`envto`) plus the hit strand:
```bash
arc.tblout.gz
bac.tblout.gz
```

3) **Run metadata**

Captures what was run, where, and when:
```bash
about_alignment.txt
```
Records at least:
- Timestamp
- Input file paths (HMMs + staged inputs)
- Output directory paths
- Tool versions (first line of `nhmmer -h` / `hmmalign -h`)
- The exact commands executed (or CLI args used)

4) **Logs**
```bash
arc.hmmer.log
bac.hmmer.log
```

**Implementation notes:**
- **No chunking:** Step 2 runs each HMM search once over the full `all_asvs.fna` (no splitting into chunks).
- **No resume logic:** Step 2 is a simple bash wrapper; it overwrites outputs in the `02_hmmer/` directory.
- **No line wrapping:** runs `nhmmer` with `--notextw` so the human-readable log output is not wrapped.
- **No alignment blocks in logs:** runs `nhmmer` with `--noali` so the logs exclude per-hit alignments (keeps log files small).
- **E-value reporting threshold:** runs `nhmmer` with `-E ${NHMMER_REPORT_MAX_EVALUE}` (and optionally `--incE ${NHMMER_REPORT_MAX_EVALUE}`) so that Step 2 does not write huge tblout files full of weak hits. This threshold is about keeping the raw output *manageable*, not about defining what counts as a "good" hit for consensus calling.
- **RF sanity check:** after `hmmalign`, verifies that each `*.sto` contains a `#=GC RF` line. If not present, fails with a clear error message.

**Running Step 2:**
```bash
bash Scripts/asvs2truncspec_02_hmmer.sh \
  --arc_hmm InputData/arc_16s.hmm \
  --bac_hmm InputData/bac_16s.hmm \
  --staged_dir asvs2truncspec_out/intermediates/01_staged_inputs \
  --out_dir asvs2truncspec_out/intermediates/02_hmmer \
  --threads 18 \
  --evalue 1e-3 \
  --verbose
```

CLI options:
- `--arc_hmm PATH`: Path to archaeal 16S HMM file (required)
- `--bac_hmm PATH`: Path to bacterial 16S HMM file (required)
- `--staged_dir PATH`: Path to Step 1 staged inputs directory (required)
- `--out_dir PATH`: Path to output directory for HMMER results (required)
- `--threads N`: Number of threads for nhmmer (default: 1)
- `--evalue E`: E-value threshold for nhmmer reporting (default: 1e-3)
- `--verbose`: Print verbose output

### Processing Step 3 (asvs2truncspec_03_agg.py)

This step converts the raw HMMER outputs into the actual `.truncspec` entries (separately for the archaeal and bacterial reference sequences).

Conceptually, Step 3 does three things:
1) Parse the `*.sto` reference alignments to create a mapping of **model position → reference position**.
2) Parse the `*.tblout` files to get the best hit (if any) per ASV to each HMM, and convert `hmmfrom/hmmto` into `ref_start/ref_end`.
3) Aggregate per dataset and compute a robust consensus region.

**Model → reference mapping (and how we handle reference gaps)**

Step 3 translates HMMER **model coordinates** (`hmmfrom/hmmto`) into **reference-sequence coordinates** (so we can write `.truncspec` entries).

A `model_pos (1-based) -> ref_pos (1-based)` mapping is created using the reference `*.sto` file:
- Walk the Stockholm alignment left-to-right, concatenating all blocks for the reference sequence and the `#=GC RF` line (ignore all other lines for mapping).
- Use the `#=GC RF` annotation (present in the provided 16S HMMs) to define which alignment columns correspond to **HMM match states** (i.e. "real" model positions). Convention: `x` = match-state column, `.` = non-match column.
- Maintain two counters:
  - `model_pos`: increments only on match-state columns
  - `ref_pos`: increments only when the reference sequence has a nucleotide (A/C/G/T/U/N; case-insensitive) in that column
- For each match-state column:
  - If the reference has a nucleotide: record `model_pos -> ref_pos`
  - If the reference has a gap at that column (`-` or `.`): record **no mapping** for that `model_pos` (leave it missing)

When converting an ASV hit's `hmmfrom/hmmto` to reference coordinates, missing mappings are handled by "snapping" to the nearest mapped model position:
- `ref_start`: start at `hmmfrom` and move **forward** until a mapped `model_pos` is found.
- `ref_end`: start at `hmmto` and move **backward** until a mapped `model_pos` is found.
- If snapping runs off the ends (no mapped position found), treat that ASV as having no usable coordinate for that domain (record as missing for QC/consensus).

This approach is intentionally "close enough" and robust: if the reference has a short gap run at a boundary, the snapping error is bounded by that local gap run length (typically well within ~5–20 bp tolerance).

**Inputs:**
- `INFO_OUT_DIR/intermediates/01_staged_inputs/dataset_manifest.tsv`
- `INFO_OUT_DIR/intermediates/02_hmmer/arc_ref.sto`
- `INFO_OUT_DIR/intermediates/02_hmmer/bac_ref.sto`
- `INFO_OUT_DIR/intermediates/02_hmmer/arc.tblout.gz`
- `INFO_OUT_DIR/intermediates/02_hmmer/bac.tblout.gz`

**Outputs:**

1) **Final `.truncspec`**

Written to `TRUNCSPEC_OUT_PATH`, for example:
```truncspec
ARC_REF_SEQ_ID: RS_GCF_022846175.1~NZ_AP025587.1-#2
BAC_REF_SEQ_ID: RS_GCF_030545895.1~NZ_JAUOMX010000042.1
AGP: arc_start=515, arc_end=716, bac_start=412, bac_end=637, min_len=152, max_len=276
#HMC_Dataset_001: arc_start=NA, arc_end=NA, bac_start=412, bac_end=637, min_len=152, max_len=276
```

2) **Intermediate mappings and QC summaries** (under `INFO_OUT_DIR/intermediates/03_results/`)

These are intended to make debugging easy without re-running HMMER:
```bash
arc_model_to_ref.tsv
bac_model_to_ref.tsv
asv_hits_arc.tsv
asv_hits_bac.tsv
dataset_region_calls.tsv
dataset_region_qc.tsv
hmmer_filter_stats.tsv
cross_domain_bootstrapping.tsv
region_renamings.tsv
```

The `cross_domain_bootstrapping.tsv` file is only written if bootstrapping is enabled. The `region_renamings.tsv` file is only written if region renaming is enabled.

**Hit filtering and selection logic:**

Since we do not know whether ASVs are archaeal or bacterial, we allow hits to both HMMs, but we apply a best-HMM gate per ASV to prevent "smeared" consensus. Archaeal coordinates are computed from archaeal hits that pass the gate, and bacterial coordinates are computed from bacterial hits that pass the gate.

**Important hit selection rule (avoid "smeared" consensus):**
- For each ASV, select **a single best domain hit** to a given HMM using **bitscore** (tblout `score`; highest wins).

**Basic hit filtering rule:**
- Before selecting the "best" hit for an ASV, filter candidate hits by:
  - `E-value <= ARC_HMMER_MAX_EVALUE` (archaea) or `E-value <= BAC_HMMER_MAX_EVALUE` (bacteria), both tighter than Step 2's `NHMMER_REPORT_MAX_EVALUE`
  - `coverage_frac >= MIN_HIT_COVERAGE_FRAC`
- Define `coverage_frac` as the fraction of the ASV sequence covered by the selected domain hit, computed from tblout coordinates:
  - `coverage_frac = (abs(envto - envfrom) + 1) / asv_length`
  - Use the tblout `sq len` column for `asv_length` (no FASTA re-parse required).
- **Strand handling:** allow both `+` and `-` hits; do not error on `-`. For QC, record strand counts per HMM. Use `min(envfrom, envto)` / `max(envfrom, envto)` when storing any hit span coordinates.

**Best-HMM gate:**
- After selecting the best hit per HMM for an ASV, if the ASV has passing hits to both HMMs, keep only the HMM with the higher bitscore (tblout `score`).
- Tie-breakers: lower `E-value` wins; if still tied, prefer bacterial for determinism.

**Deterministic consensus rule:**
- For each dataset and each HMM (arc/bac), compute a consensus `start` and `end` from the filtered per-ASV hit coordinates.
- Use the median (middle value after sorting; if even, take the lower of the two middle values) to get `start_med` and `end_med`.
- Drop any ASV hit where `abs(start - start_med) > OUTLIER_TOL_BP` **or** `abs(end - end_med) > OUTLIER_TOL_BP`.
- Recompute the medians on the remaining hits and use those as the final consensus `start` and `end`.
- Require at least `ARC_MIN_HITS_PER_DATASET` (arc) or `BAC_MIN_HITS_PER_DATASET` (bac) hits **after** outlier filtering; otherwise mark that domain as `NA` for that dataset.

Datasets with too few supporting hits (or with highly inconsistent hits) get `NA` coordinates for that domain **in the QC tables**, and the corresponding `.truncspec` line is commented out so it is not used by `extract16s.sh`.

**Optional cross-domain bootstrapping:**
- Purpose: fill missing domain coords when a dataset has only one domain present.
- Reference set: only datasets that already have **both** domains (arc and bac) are eligible as reference regions.
- For a missing domain coordinate (start or end), use the known domain coordinate to locate the closest reference region by that same coordinate:
  - Example: if arc is missing and bac_start is known, find the reference region with bac_start closest to this bac_start.
  - If the absolute distance is `> CROSS_DOMAIN_BOOTSTRAPPING_MAX_DISTANCE`, fail that coordinate and log it.
  - If the distance is within the max, bootstrap using: `bootstrapped_coord = ref_missing_coord + distance`.
- Start and end are bootstrapped independently (they may use different reference regions).
- If there is only one dataset total, bootstrapping is disabled with a warning.
- After bootstrapping, recompute `min_len` / `max_len` using the same rules as usual.

**Optional region redundancy minimisation:**
- Purpose: merge near-identical regions across datasets to reduce redundancy.
- Only datasets with complete coordinates in both domains are eligible for merging.
- Two regions are similar if **all four** coordinate differences are within `REGION_REDUNDANCY_MINIMISATION_MAX_DISTANCE`:
  - `abs(arc_start_a - arc_start_b)`
  - `abs(arc_end_a - arc_end_b)`
  - `abs(bac_start_a - bac_start_b)`
  - `abs(bac_end_a - bac_end_b)`
- Merging uses a simple greedy grouping that requires a candidate to be within distance of every member of a group.
- Merged region names are formed by joining dataset IDs with `REGION_REDUNDANCY_MINIMISATION_NAME_JOINER` (e.g. `HMC_Dataset_001+HMC_Dataset_014`).
- New coordinates are the mean of merged coordinates, rounded to the nearest integer.
- After merging, recompute `min_len` / `max_len` using the same rules as usual.

**Optional region renaming:**
- Purpose: rename regions to simpler, more uniform names (especially useful after redundancy minimisation creates long concatenated names).
- Three renaming methods are available:
  - `uniform`: Rename all regions to `Reg-001`, `Reg-002`, etc. using a global counter.
  - `dataset`: Rename regions to `{dataset}-001`, `{dataset}-002`, etc. where `{dataset}` is extracted from the region name. Uses a global counter.
  - `region`: Rename regions to `{region}-001`, `{region}-002`, etc. where `{region}` is the 16S region (e.g. "V4", "V3-V4") from the `region` column in `dataset_manifest.tsv`. Each unique region has its own counter starting from 001 (e.g. "V4-001", "V4-002", "V3-V4-001", "V3-V4-002").
- For the `dataset` method:
  - If the region name contains `REGION_REDUNDANCY_MINIMISATION_NAME_JOINER` (e.g. `+`), the name is split and the first dataset name is used as the prefix.
  - The global counter ensures that all region names remain unique across all datasets.
- For the `region` method:
  - The 16S region is determined from the `region` column in `dataset_manifest.tsv`.
  - If region redundancy minimisation is enabled, all datasets in a merged region must have the same 16S region; otherwise an error is raised.
  - Requires that all regions have valid region metadata (not "UNK"); otherwise an error is raised.
- The renaming produces:
  - A new `.truncspec` file with renamed regions (the final output).
  - A `_no-renaming.truncspec` file preserving the original names.
  - A `region_renamings.tsv` file mapping old names to new names.
- Renamed regions are used in the final plots.

Finally, this step calculates `min_len` and `max_len` using the existing logic (using both domains if available, otherwise just one):
- `min_len = min(found_lengths) - MIN_LEN_BUFFER`
- `max_len = max(found_lengths) + MAX_LEN_BUFFER`

**Cross-domain tuning stats:**
- Writes a summary table (`hmmer_filter_stats.tsv`) with per-HMM counts and rates:
  - total ASVs searched
  - total hits reported
  - hits passing `E-value` + coverage filters
  - hits retained after best-HMM gate
  - percent of ASVs with a passing hit (pre- and post-gate)
  - strand counts (`+`/`-`)
- These stats make it easy to tune HMMER filtering so that archaeal hit rates are low (e.g. 1-5%) and bacterial hit rates are high (e.g. ~85%), without changing the core logic.

**Running Step 3:**
```bash
python Scripts/asvs2truncspec_03_agg.py
```

### Resulting Outputs

After running all steps, the main outputs are:
- `TRUNCSPEC_OUT_PATH` (the `.truncspec` file used by `extract16s.sh`)
- `INFO_OUT_DIR/logs/` (high-level run logs)
- `INFO_OUT_DIR/intermediates/` (staged inputs, raw HMMER outputs, and QC tables)

The intended workflow is that if you adjust input ASVs or dataset metadata, you re-run Step 1 → Step 2 → Step 3. If you are only tweaking the consensus/QC logic, you can usually re-run Step 3 only.


---


## Usage Guide

### Running the Complete Pipeline

To run the complete `asvs2truncspec` pipeline:

1. **Configure the scripts:**
   - Edit `Scripts/asvs2truncspec_01_prep.py` to set your paths and reference IDs
   - Edit `Scripts/asvs2truncspec_03_agg.py` to set your filtering thresholds and optional features

2. **Run Step 1 (Staging):**
   ```bash
   python Scripts/asvs2truncspec_01_prep.py
   ```

3. **Run Step 2 (HMMER searches):**
   ```bash
   bash Scripts/asvs2truncspec_02_hmmer.sh \
     --arc_hmm InputData/arc_16s.hmm \
     --bac_hmm InputData/bac_16s.hmm \
     --staged_dir asvs2truncspec_out/intermediates/01_staged_inputs \
     --out_dir asvs2truncspec_out/intermediates/02_hmmer \
     --threads 18 \
     --evalue 1e-3 \
     --verbose
   ```

4. **Run Step 3 (Process outputs):**
   ```bash
   python Scripts/asvs2truncspec_03_agg.py
   ```

5. **Review outputs:**
   - Check the final `.truncspec` file at `TRUNCSPEC_OUT_PATH`
   - Review QC tables in `INFO_OUT_DIR/intermediates/03_results/`
   - Check `hmmer_filter_stats.tsv` to ensure appropriate hit rates

### Iterating on Parameters

If you need to adjust filtering or consensus parameters:
- Edit the configuration constants in `Scripts/asvs2truncspec_03_agg.py`
- Re-run Step 3 only (no need to re-run Steps 1 and 2)

If you need to adjust HMMER E-value thresholds:
- Edit `NHMMER_REPORT_MAX_EVALUE` in Step 2
- Re-run Steps 2 and 3

### Troubleshooting

**Few or no archaeal hits:**
- Archaeal ASVs are typically much less common than bacterial
- Check `hmmer_filter_stats.tsv` to see hit rates
- Consider adjusting `ARC_HMMER_MAX_EVALUE` or `ARC_MIN_HITS_PER_DATASET`

**Many incomplete regions:**
- Check `dataset_region_qc.tsv` to see which datasets failed and why
- Enable cross-domain bootstrapping if one domain consistently has coordinates
- Adjust `MIN_HITS_PER_DATASET` if you have very small datasets

**Regions appear too similar:**
- Enable region redundancy minimisation
- Adjust `REGION_REDUNDANCY_MINIMISATION_MAX_DISTANCE` to control merge threshold


---


## Miscellaneous Details

- Typical ASV length range across datasets is between 90-500 bp
- We never expect mixed regions within the same dataset (beyond rare outliers)
- Cross-domain bootstrapping requires at least 2 datasets to function
