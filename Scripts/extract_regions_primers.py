"""
extract_regions_primers.py 

A demo script for extracting variable regions from 16S rRNA sequences **using primers**

Overview:
 - Filters sequences by length and ambiguous base content.
 - If REQUIRE_ALL_REGIONS is True, only keep sequences that pass through all regions.
 - This script should be used on the output of the filter_database.py script.

Input: 16S GTDB database file (e.g. ssu_all_r220_filtered.fna)

Output: Directory with results:
    - [NAME]_seqs.fasta: Extracted sequences
    - about.txt: Summary of processing details
    - FULL_seqs.fasta: Full length sequences

Usage
 - This is a standalone python script
 - You just need to input the primers (PRIMER_PAIRS) and the database file (INPUT_FILE) plus any other parameters
 - The script will output the extracted sequences into the output directory
"""

# Imports
import re
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
import itertools
import statistics
import os
import shutil
import datetime
import matplotlib.pyplot as plt

# Input database file with full length sequences
INPUT_FILE = "D:/16S_databases/filtered_databases/ssu_all_r220_filtered.fna"

# Output
OUTPUT_DIR = "D:/16S_databases/micro16S_databases/test_db/"
OUTPUT_DB_PATH = OUTPUT_DIR + "[NAME]_seqs.fasta"
ABOUT_PATH = OUTPUT_DIR + "about.txt"
SEQ_LENGTHS_PLOTS = OUTPUT_DIR + "plots/[NAME]_seqs_lengths.png"

# File to copy to the output directory (unless None)
TAXON_NAME_TO_DISTANCE_BETWEEN_ANCESTORS_DICT_FILE_IN = "D:/16S_databases/GTDBtk/taxon_to_dist_dict.pkl"
TAXON_NAME_TO_DISTANCE_BETWEEN_ANCESTORS_DICT_FILE_OUT = OUTPUT_DIR + "taxon_to_dist_dict.pkl"


# Number of CPU cores to use
NUM_PROCESSORS = 18

# Primer pairs to process (each results in a separate output file)
PRIMER_PAIRS = {
    "V4": {
        "fwd_primers": [
            {'name': '515F (Parada)', 'seq': 'GTGYCAGCMGCCGCGGTAA'},
            # {'name': '515F (Caporaso)', 'seq': 'GTGCCAGCMGCCGCGGTAA'},
        ],
        "rev_primers": [
            {'name': '806R (Apprill)', 'seq': 'GGACTACNVGGGTWTCTAAT'},
            # {'name': '806R (Caporaso)', 'seq': 'GGACTACHVGGGTWTCTAAT'},
        ],
        "min_length": 245,
        "max_length": 260,
        "max_bad_chars": 3,
        "sources": [
            "https://earthmicrobiome.org/protocols-and-standards/16s/",
            "http://doi.org/10.1111/1462-2920.13023",
            "http://doi.org/10.3354/ame01753"
        ]
    },
    "V3-V4": {
        "fwd_primers": [
            {'name': 'Bakt_341F', 'seq': 'CCTACGGGNGGCWGCAG'},
        ],
        "rev_primers": [
            {'name': 'Bakt_805R', 'seq': 'GACTACHVGGGTATCTAATCC'},
        ],
        "min_length": 385,
        "max_length": 445,
        "max_bad_chars": 3,
        "sources": [
            "https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2019.01937/full"
        ]
    }
}

# If True, only keep sequences that pass filters for all regions
REQUIRE_ALL_REGIONS = True

# If True, write the sequence indices to the output file. For example: 
    # >{0}RS_GCF_008416195.1~NZ_VVXS01000056.1 d__Bacteria;p__Bacteroidota;c__Bacteroidia;...
    # >{1}RS_GCF_015238635.1~NZ_CP049958.1-#7 d__Bacteria;p__Bacillota_A;c__Clostridia;...
WRITE_SEQ_INDICES = True

# Can't write sequence indices if not requiring all regions
if WRITE_SEQ_INDICES and not REQUIRE_ALL_REGIONS:
    raise ValueError("Can't write sequence indices if not requiring all regions")

# IUPAC ambiguity codes lookup
AMBITABLE = {
    'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C',
    'R': '[AG]', 'Y': '[CTU]', 'K': '[GTU]', 'M': '[AC]',
    'S': '[CG]', 'W': '[ATU]', 'B': '[CGTU]', 'D': '[AGTU]',
    'H': '[ACTU]', 'V': '[ACG]', 'N': '[ACGTU]', 'X': 'X',
    '-': '-'
}
# Add lowercase versions
AMBITABLE.update({k.lower(): v for k, v in AMBITABLE.items()})

def get_parent_directory(path):
    return os.path.dirname(path.rstrip("/").rstrip("\\"))

def expand_sequence(seq):
    """Converts a sequence with IUPAC codes to regex pattern"""
    return "".join(AMBITABLE.get(base, base) for base in seq)

def load_primers(fwd_primers, rev_primers):
    """Converts lists of primer entries into combined regex patterns for forward and reverse primers"""

    # Build combined forward primer pattern
    fwd_patterns = []
    for primer in fwd_primers:
        seq_value = primer['seq']
        variations = (
            [seq_value] +
            [f"{seq_value[:x]}N{seq_value[x:]}" for x in range(1, len(seq_value))] +  # insertions
            [f"{seq_value[:x]}N{seq_value[x+1:]}" for x in range(len(seq_value)-1)] +  # substitutions
            [f"{seq_value[:x]}{seq_value[x+1:]}" for x in range(len(seq_value))]      # deletions
        )
        fwd_patterns.extend(expand_sequence(v) for v in variations)
    combined_fwd_pattern = "(" + "|".join(fwd_patterns) + ")"

    # Build combined reverse primer pattern (including reverse complement variants)
    rev_patterns = []
    for primer in rev_primers:
        seq_value = primer['seq']
        for base in [seq_value, str(Seq(seq_value).reverse_complement())]:
            variations = (
                [base] +
                [f"{base[:x]}N{base[x:]}" for x in range(1, len(base))] +
                [f"{base[:x]}N{base[x+1:]}" for x in range(len(base)-1)] +
                [f"{base[:x]}{base[x+1:]}" for x in range(len(base))]
            )
            rev_patterns.extend(expand_sequence(v) for v in variations)
    combined_rev_pattern = "(" + "|".join(rev_patterns) + ")"

    return (combined_fwd_pattern, combined_rev_pattern)

def process_record(args):
    """Process a single sequence record, applying filters and trimming"""
    record, primers, min_len, max_len, badchars = args
    
    # Filter sequences with too many ambiguous bases
    seq = str(record.seq)
    if badchars > -1 and len(re.findall(r'[rykmswbdhvnx-]', seq, re.IGNORECASE)) > badchars:
        return None
        
    # Find and trim between primers
    fp = re.search(primers[0], seq, re.IGNORECASE)
    if not fp:
        return None
    seq = seq[fp.span()[1]:]
        
    rp = re.search(primers[1], seq, re.IGNORECASE)
    if not rp:
        return None
    seq = seq[:rp.span()[0]]
    
    if not min_len <= len(seq) <= max_len:
        return None

    record.seq = Seq(seq)
    return record

def process_file(input_file, num_processes, primers, min_length, max_length, max_bad_chars):
    """Process sequences in parallel using multiprocessing"""
    sequence_lengths = []
    total_count = 0
    failed_count = 0
    processed_records = []  # Store successful records
    
    with Pool(num_processes) as pool:
        records = ((rec, primers, min_length, max_length, max_bad_chars)
                   for rec in SeqIO.parse(input_file, "fasta"))
        while True:
            chunk = list(itertools.islice(records, num_processes))
            if not chunk:
                break
            total_count += len(chunk)
            results = pool.map(process_record, chunk)
            # Count successful results and track lengths
            for res in results:
                if res:
                    sequence_lengths.append(len(res.seq))
                    processed_records.append(res)
                else:
                    failed_count += 1
    
    # Calculate statistics
    if sequence_lengths:
        stats = {
            'total': total_count,
            'failed': failed_count,
            'min_len': min(sequence_lengths),
            'p1': statistics.quantiles(sequence_lengths, n=100)[0],
            'p10': statistics.quantiles(sequence_lengths, n=10)[0],
            'median': statistics.median(sequence_lengths),
            'p90': statistics.quantiles(sequence_lengths, n=10)[-1],
            'p99': statistics.quantiles(sequence_lengths, n=100)[-1],
            'max_len': max(sequence_lengths)
        }
    else:
        stats = {
            'total': total_count,
            'failed': failed_count,
            'min_len': 0, 'p1': 0, 'p10': 0, 'median': 0, 'p90': 0, 'p99': 0, 'max_len': 0
        }
    return stats, processed_records

def plot_sequence_lengths(lengths_dict, output_path):
    """Create histogram plots of sequence lengths for each region.
    
    Args:
        lengths_dict: Dictionary with region names as keys and lists of lengths as values
        output_path: Path template for saving plots (should contain [NAME] placeholder)
    """
    # Create plots directory if it doesn't exist
    os.makedirs(get_parent_directory(output_path), exist_ok=True)
    
    # Plot each region
    for region_name, lengths in lengths_dict.items():
        plt.figure(figsize=(10, 6))
        plt.hist(lengths, bins=50, edgecolor='black')
        plt.title(f'Sequence Length Distribution - {region_name} Region')
        plt.xlabel('Sequence Length (bp)')
        plt.ylabel('Count')
        
        # Add vertical lines for min/max length filters from PRIMER_PAIRS
        min_len = PRIMER_PAIRS[region_name]['min_length']
        max_len = PRIMER_PAIRS[region_name]['max_length']
        plt.axvline(x=min_len, color='r', linestyle='--', label=f'Min Length ({min_len} bp)')
        plt.axvline(x=max_len, color='r', linestyle='--', label=f'Max Length ({max_len} bp)')
        
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # Save plot
        output_file = output_path.replace("[NAME]", region_name)
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

if __name__ == '__main__':
    # Store all stats and records for each region
    all_stats = {}
    all_records = {}
    full_records = []
    sequence_lengths_by_region = {}

    # If the directory doesn't exist, create it
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load full sequences first
    print("\nLoading full sequences...")
    full_records = list(SeqIO.parse(INPUT_FILE, "fasta"))

    # Get the number of sequences in the full database
    full_seq_count = len(full_records)
    
    # Process each primer pair
    for region_name, primer_pair in PRIMER_PAIRS.items():

        print(f"\nProcessing {region_name} region...")
        print("Forward primers:")
        for p in primer_pair['fwd_primers']:
            print(f"  {p['name']}: {p['seq']}")
        print("Reverse primers:")
        for p in primer_pair['rev_primers']:
            print(f"  {p['name']}: {p['seq']}")

        # Load the primers (combined forward and reverse patterns)
        primers = load_primers(primer_pair['fwd_primers'], primer_pair['rev_primers'])
        
        # Unpack the dictionary into separate variables
        min_length = primer_pair['min_length']
        max_length = primer_pair['max_length']
        max_bad_chars = primer_pair['max_bad_chars']
        
        # Process sequences but don't write output yet
        with open(INPUT_FILE) as input_file:
            stats, records = process_file(input_file, NUM_PROCESSORS, primers, min_length, max_length, max_bad_chars)
        
        # Store stats and records for this region
        all_stats[region_name] = {
            'stats': stats,
            'primers': primer_pair
        }
        all_records[region_name] = {rec.id: rec for rec in records}
        
        # Store sequence lengths for plotting
        sequence_lengths_by_region[region_name] = [len(rec.seq) for rec in records]
        
        print(f"Processed {region_name}")
        print(f"Failed sequences: {stats['failed']} out of {stats['total']} ({stats['failed']/stats['total']*100:.1f}%)")
        print("Sequence lengths:")
        print(f"  Min: {stats['min_len']}")
        print(f"  1st percentile: {stats['p1']}")
        print(f"  10th percentile: {stats['p10']}")
        print(f"  Median: {stats['median']}")
        print(f"  90th percentile: {stats['p90']}")
        print(f"  99th percentile: {stats['p99']}")
        print(f"  Max: {stats['max_len']}")

    # Create length distribution plots
    plot_sequence_lengths(sequence_lengths_by_region, SEQ_LENGTHS_PLOTS)

    # Filter sequences if REQUIRE_ALL_REGIONS is True
    if REQUIRE_ALL_REGIONS:
        # Find sequences present in all regions
        common_ids = set.intersection(*[set(records.keys()) for records in all_records.values()])
        
        # Get the order of sequences from the full database to maintain consistency
        ordered_ids = [rec.id for rec in full_records if rec.id in common_ids]
        
        # Write filtered outputs in consistent order
        for region_name, records in all_records.items():
            output_file = OUTPUT_DB_PATH.replace("[NAME]", region_name)
            with open(output_file, "w") as out_file:
                i = 0
                for _, seq_id in enumerate(ordered_ids):
                    record = records[seq_id]
                    index_str = "{" + str(i) + "}" if WRITE_SEQ_INDICES else ""
                    out_file.write(f">{index_str}{record.description}\n{str(record.seq)}\n")
                    i += 1
        
        # Write filtered full sequences in same order
        full_output = OUTPUT_DB_PATH.replace("[NAME]", "FULL")
        with open(full_output, "w") as out_file:
            i = 0
            for _, rec in enumerate(full_records):
                if rec.id in common_ids:
                    index_str = "{" + str(i) + "}" if WRITE_SEQ_INDICES else ""
                    out_file.write(f">{index_str}{rec.description}\n{str(rec.seq)}\n")
                    i += 1
        
        # Calculate and display statistics
        final_count = len(common_ids)
        percent_remaining = (final_count / full_seq_count) * 100
        print("\nAfter requiring sequences to be present in all regions:")
        print(f"Final sequence count: {final_count:,} ({percent_remaining:.1f}% of original)")
    else:
        # Write unfiltered outputs
        for region_name, records in all_records.items():
            output_file = OUTPUT_DB_PATH.replace("[NAME]", region_name)
            with open(output_file, "w") as out_file:
                for record in records.values():
                    out_file.write(f">{record.description}\n{str(record.seq)}\n")
        # Copy full sequences
        shutil.copy(INPUT_FILE, OUTPUT_DB_PATH.replace("[NAME]", "FULL"))

    # Write the about file
    with open(ABOUT_PATH, 'w') as f:
        f.write("16S Variable Region Extraction Summary\n")
        f.write("====================================\n\n")
        f.write(f"Input database: {INPUT_FILE}\n")
        f.write(f"Processing date: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

        for region_name, data in all_stats.items():
            header = f"{region_name} Region"
            f.write(header + "\n")
            f.write("-" * len(header) + "\n")
            
            # Primer information for multiple primers
            f.write("Forward primers:\n")
            for p in data['primers']['fwd_primers']:
                f.write(f"  {p['name']}: {p['seq']}\n")
            f.write("Reverse primers:\n")
            for p in data['primers']['rev_primers']:
                f.write(f"  {p['name']}: {p['seq']}\n")
            f.write(f"Length filters: {data['primers']['min_length']}-{data['primers']['max_length']} bp\n")
            f.write(f"Max ambiguous bases: {data['primers']['max_bad_chars']}\n\n")
            
            # Sources
            f.write("Sources:\n")
            for source in data['primers']['sources']:
                f.write(f"- {source}\n")
            f.write("\n")
            
            # Statistics
            stats = data['stats']
            f.write("Results:\n")
            f.write(f"Total sequences processed: {stats['total']}\n")
            f.write(f"Failed sequences: {stats['failed']} ({stats['failed']/stats['total']*100:.1f}%)\n")
            f.write("Sequence length distribution:\n")
            f.write(f"  Minimum: {stats['min_len']} bp\n")
            f.write(f"  1st percentile: {stats['p1']} bp\n")
            f.write(f"  10th percentile: {stats['p10']} bp\n")
            f.write(f"  Median: {stats['median']} bp\n")
            f.write(f"  90th percentile: {stats['p90']} bp\n")
            f.write(f"  99th percentile: {stats['p99']} bp\n")
            f.write(f"  Maximum: {stats['max_len']} bp\n\n")

        if REQUIRE_ALL_REGIONS:
            # Add statistics to about file
            f.write("\nSequences Present in All Regions\n")
            f.write("==============================\n")
            f.write(f"Final sequence count: {final_count:,}\n")
            f.write(f"Percentage of original database: {percent_remaining:.1f}%\n")
    
    # Copy the taxon name to distance file if specified
    if TAXON_NAME_TO_DISTANCE_BETWEEN_ANCESTORS_DICT_FILE_IN:
        shutil.copy(TAXON_NAME_TO_DISTANCE_BETWEEN_ANCESTORS_DICT_FILE_IN, TAXON_NAME_TO_DISTANCE_BETWEEN_ANCESTORS_DICT_FILE_OUT)

    print("Done!")


