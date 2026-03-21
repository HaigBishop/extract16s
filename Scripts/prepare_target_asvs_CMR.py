#!/usr/bin/env python3
"""
This script is specifically for the Celiac Microbiome Repository (CMR). See more information here: https://github.com/CeliacMicrobiomeRepo/celiac-repository

Script to prepare ASV files for asvs2truncspec.
It copies seqs.fna files from the celiac-repository to the target_asvs directory
and generates the corresponding _datasets.tsv metadata files.
"""

import os
import csv

# Configuration Paths
CELIAC_REPO_DIR = os.path.expanduser("~/Repos/celiac-repository")
EXTRACT16S_REPO_DIR = os.path.expanduser("~/Repos/extract16s")

SOURCE_DATASETS_DIR = os.path.join(CELIAC_REPO_DIR, "16S_datasets")
INCLUDED_DATASETS_TSV = os.path.join(CELIAC_REPO_DIR, "included_datasets.tsv")
TARGET_ASVS_DIR = os.path.join(EXTRACT16S_REPO_DIR, "target_asvs")

def load_dataset_regions(tsv_path):
    """
    Loads dataset regions from the included_datasets.tsv file.
    Returns a dictionary mapping Dataset_ID to Amplicon_Region.
    """
    dataset_regions = {}
    print(f"Reading dataset metadata from: {tsv_path}")
    
    try:
        with open(tsv_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                dataset_id = row.get('Dataset_ID')
                region = row.get('Amplicon_Region')
                
                if dataset_id:
                    # Clean up region string if necessary
                    if not region or region.strip() == '':
                        region = 'UNK'
                    else:
                        region = region.strip()
                        
                    dataset_regions[dataset_id] = region
    except Exception as e:
        print(f"Error reading {tsv_path}: {e}")
        return {}
        
    print(f"Loaded regions for {len(dataset_regions)} datasets.")
    return dataset_regions

def main():
    # Ensure target directory exists
    if not os.path.exists(TARGET_ASVS_DIR):
        print(f"Creating target directory: {TARGET_ASVS_DIR}")
        os.makedirs(TARGET_ASVS_DIR, exist_ok=True)
    else:
        print(f"Target directory exists: {TARGET_ASVS_DIR}")

    # Load region metadata
    dataset_regions = load_dataset_regions(INCLUDED_DATASETS_TSV)
    
    if not dataset_regions:
        print("Warning: No dataset regions loaded. Proceeding with 'UNK' for unknown regions.")

    # Iterate through datasets in the source directory
    if not os.path.exists(SOURCE_DATASETS_DIR):
        print(f"Error: Source directory not found: {SOURCE_DATASETS_DIR}")
        return

    print(f"Scanning for datasets in: {SOURCE_DATASETS_DIR}")
    
    processed_count = 0
    skipped_count = 0
    
    for item in sorted(os.listdir(SOURCE_DATASETS_DIR)):
        dataset_dir = os.path.join(SOURCE_DATASETS_DIR, item)
        
        # Check if it's a directory and likely a dataset (starts with 16S_)
        if os.path.isdir(dataset_dir) and item.startswith("16S_"):
            dataset_id = item
            src_fna = os.path.join(dataset_dir, "seqs.fna")
            
            if os.path.exists(src_fna):
                # Determine region
                region = dataset_regions.get(dataset_id)
                
                if region is None:
                    print(f"Warning: Region not found for {dataset_id} in metadata. Using 'UNK'.")
                    region = "UNK"
                
                # Define target paths
                dest_fna = os.path.join(TARGET_ASVS_DIR, f"{dataset_id}.fna")
                dest_tsv = os.path.join(TARGET_ASVS_DIR, f"{dataset_id}_datasets.tsv")
                
                # Parse and rename ASVs, then write to target FASTA and TSV
                try:
                    asv_ids = []
                    unique_asv_ids = []
                    
                    with open(src_fna, 'r', encoding='utf-8') as f_in, \
                         open(dest_fna, 'w', encoding='utf-8') as f_out:
                        for line in f_in:
                            if line.startswith('>'):
                                # Extract original ID
                                original_header = line[1:].strip()
                                original_id = original_header.split()[0]
                                
                                # Create new unique ID
                                unique_id = f"{dataset_id}_{original_id}"
                                
                                # Write renamed header
                                # Preserve description if any existed
                                parts = original_header.split(None, 1)
                                if len(parts) > 1:
                                    f_out.write(f">{unique_id} {parts[1]}\n")
                                else:
                                    f_out.write(f">{unique_id}\n")
                                
                                unique_asv_ids.append(unique_id)
                            else:
                                f_out.write(line)
                    
                    if not unique_asv_ids:
                        print(f"Warning: No ASVs found in {src_fna}")
                    
                    # Create _datasets.tsv file using unique IDs
                    with open(dest_tsv, 'w', encoding='utf-8', newline='') as f:
                        writer = csv.writer(f, delimiter='\t')
                        writer.writerow(['ASV_ID', 'Dataset_ID', 'Region'])
                        for unique_id in unique_asv_ids:
                            writer.writerow([unique_id, dataset_id, region])
                            
                except Exception as e:
                    print(f"Error processing {dataset_id}: {e}")
                    continue
                
                processed_count += 1
                if processed_count % 10 == 0:
                    print(f"Processed {processed_count} datasets...")
            else:
                # print(f"Skipping {dataset_id}: seqs.fna not found.")
                skipped_count += 1
    
    print("-" * 30)
    print(f"Processing complete.")
    print(f"Total datasets processed: {processed_count}")
    print(f"Datasets skipped (no seqs.fna): {skipped_count}")
    print(f"Output directory: {TARGET_ASVS_DIR}")

if __name__ == "__main__":
    main()
