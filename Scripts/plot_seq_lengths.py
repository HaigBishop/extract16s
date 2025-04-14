"""
Generates sequence length histogram plots for FASTA files.

Reads FASTA files matching the pattern {TARGET_DIR}/*_seqs.fasta.
For each input file, it calculates the lengths of all sequences and
generates a histogram plot.

The plots are saved as PNG files in the {TARGET_DIR}/seq_len_plots/
directory, named like *_seqs_lengths.png. The output directory
is created if it does not exist.

Configuration is done via constants at the top of the script.
"""


# --- Imports ---
import os
import glob
from Bio import SeqIO
import matplotlib.pyplot as plt



# --- Constants ---
TARGET_DIR = "./Output/"
SHOW_PLOTS = False # Display plots interactively?
# Histogram bin widths for sequence lengths
BIN_SIZE_FULL_SEQ = 20       # For files with "FULL" in filename
BIN_SIZE_OTHER_SEQS = 2      # For files without "FULL" in filename



# --- Script ---
def plot_sequence_lengths(fasta_file, output_dir, bin_size, show_plot):
    """Reads a FASTA file, calculates sequence lengths, and plots a histogram."""
    lengths = []
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            lengths.append(len(record.seq))
    except FileNotFoundError:
        print(f"Error: Input file not found: {fasta_file}")
        return
    except Exception as e:
        print(f"Error reading {fasta_file}: {e}")
        return

    if not lengths:
        print(f"Warning: No sequences found in {fasta_file}")
        return

    num_sequences = len(lengths)
    base_filename = os.path.basename(fasta_file)
    output_filename = os.path.splitext(base_filename)[0] + "_lengths.png"
    output_path = os.path.join(output_dir, output_filename)

    # Determine bins
    min_len = min(lengths) if lengths else 0
    max_len = max(lengths) if lengths else bin_size
    bins = range(min_len, max_len + bin_size, bin_size)

    # Create plot
    plt.figure(figsize=(10, 6))
    plt.hist(lengths, bins=bins, edgecolor='black')
    plt.title(f'Sequence Length Distribution ({num_sequences} seqs)\n{base_filename}')
    plt.xlabel('Sequence Length (bp)')
    plt.ylabel('Number of Sequences')
    plt.grid(axis='y', alpha=0.75)

    # Save plot
    try:
        plt.savefig(output_path)
        print(f"Saved plot: {output_path}")
    except Exception as e:
        print(f"Error saving plot {output_path}: {e}")


    # Show plot if requested
    if show_plot:
        plt.show()

    plt.close() # Close the figure to free memory



if __name__ == "__main__":
    output_dir = os.path.join(TARGET_DIR, "seq_len_plots")
    input_pattern = os.path.join(TARGET_DIR, "*_seqs.fasta")

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Find and process FASTA files
    fasta_files = glob.glob(input_pattern)

    if not fasta_files:
        print(f"No FASTA files found matching pattern: {input_pattern}")
    else:
        print(f"Found {len(fasta_files)} FASTA files to process.")

        for fasta_file in fasta_files:
            bin_size = BIN_SIZE_FULL_SEQ if "FULL" in fasta_file else BIN_SIZE_OTHER_SEQS
            print(f"Processing: {fasta_file}")
            plot_sequence_lengths(fasta_file, output_dir, bin_size, SHOW_PLOTS)

        print("Finished processing all files.")
