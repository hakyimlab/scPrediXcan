#!/usr/bin/env python3
"""
Simple ctPred predictor - H5 files with gene keys
"""

import torch
import pandas as pd
import h5py
import numpy as np
import os
import sys
import time
from pathlib import Path
import argparse

def print_progress(current, total, start_time, prefix=""):
    """Simple progress display without external libraries"""
    elapsed = time.time() - start_time
    percent = 100 * (current / total) if total > 0 else 0
    bar_length = 40
    filled_length = int(bar_length * current // total)
    bar = '█' * filled_length + '░' * (bar_length - filled_length)
    
    sys.stdout.write(f'\r{prefix} |{bar}| {current}/{total} ({percent:.1f}%) {elapsed:.1f}s')
    sys.stdout.flush()
    
    if current == total:
        print()

def main():
    from ctPred_utils import ctPred
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Predict gene expression from H5 files')
    parser.add_argument('--model_path', required=True, help='Path to trained model (.pt file)')
    parser.add_argument('--input_folder', required=True, help='Folder containing .h5 files')
    parser.add_argument('--output_path', required=True, help='Output CSV file path')
    parser.add_argument('--ref_csv', required=True, help='Reference CSV file with gene mappings')
    parser.add_argument('--device', default='cpu', help='Device to use (cpu or cuda)')
    args = parser.parse_args()
    
    # Load model
    device = torch.device(args.device)
    model = ctPred().to(device)
    model.load_state_dict(torch.load(args.model_path, map_location=device))
    model.eval()
    
    # Load reference CSV
    ref_df = pd.read_csv(args.ref_csv)
    
    # Check which column contains the TSS_enformer_input format
    tss_col = None
    possible_tss_cols = ['TSS_enformer_input', 'tss_enformer_input', 'region', 'tss_region', 'chr_tss']
    
    for col in possible_tss_cols:
        if col in ref_df.columns:
            tss_col = col
            break
    
    if tss_col is None:
        for col in ref_df.columns:
            if any('chr' in str(val) for val in ref_df[col].head(10)):
                tss_col = col
                break
    
    if tss_col is None:
        print(f"Error: Could not find TSS region column in reference CSV.")
        sys.exit(1)
    
    # Create mapping from region to ensembl_gene_id
    ref_df[tss_col] = ref_df[tss_col].astype(str).str.replace('"', '').str.strip()
    
    # Make sure we have ensembl_gene_id column
    if 'ensembl_gene_id' not in ref_df.columns:
        gene_col = None
        for col in ['ensembl_gene_id', 'gene_id', 'gene', 'ENSEMBL']:
            if col in ref_df.columns:
                gene_col = col
                break
        
        if gene_col:
            ref_df['ensembl_gene_id'] = ref_df[gene_col].astype(str).str.replace('"', '').str.strip()
        else:
            print(f"Error: No gene ID column found in reference CSV.")
            sys.exit(1)
    else:
        ref_df['ensembl_gene_id'] = ref_df['ensembl_gene_id'].astype(str).str.replace('"', '').str.strip()
    
    # Create mapping dictionary
    region_to_gene = dict(zip(ref_df[tss_col], ref_df['ensembl_gene_id']))
    
    # Get all H5 files
    h5_files = list(Path(args.input_folder).glob("*.h5"))
    
    if not h5_files:
        print(f"Error: No .h5 files found in {args.input_folder}")
        sys.exit(1)
    
    # Process each file with simple progress display
    results = {}
    total_genes = 0
    start_time = time.time()
    
    print(f"Processing {len(h5_files)} samples...")
    
    for i, h5_file in enumerate(h5_files, 1):
        print_progress(i, len(h5_files), start_time, "Samples")
        
        with h5py.File(h5_file, 'r') as f:
            preds = {}
            gene_keys = list(f.keys())
            
            # Simple counter for genes in this file
            for j, gene_id in enumerate(gene_keys, 1):
                data = f[gene_id][:]
                features = np.mean(data, axis=0)
                
                X = torch.tensor(features, dtype=torch.float32).to(device)
                with torch.no_grad():
                    pred = model(X.unsqueeze(0)).item()
                preds[gene_id] = pred
            
            results[h5_file.stem] = preds
            total_genes += len(preds)
    
    print(f"\nProcessed {total_genes} gene predictions")
    
    # Create initial DataFrame
    df = pd.DataFrame.from_dict(results, orient='index').T
    
    # Process row names: remove '_predictions' suffix
    df.index = df.index.str.replace('_predictions', '')
    
    # Map regions to gene IDs and filter
    valid_indices = []
    gene_ids = []
    unmatched_count = 0
    
    total_regions = len(df.index)
    
    for i, region in enumerate(df.index, 1):
        if i % 100 == 0 or i == total_regions:
            print_progress(i, total_regions, start_time, "Final file framing")
        
        if region in region_to_gene:
            valid_indices.append(region)
            gene_ids.append(region_to_gene[region])
        else:
            unmatched_count += 1
    
    # Filter DataFrame to only include matched regions
    if valid_indices:
        df_matched = df.loc[valid_indices].copy()
        df_matched.index = gene_ids
        df_matched.index.name = 'NAME'
        df_matched = df_matched.reset_index()
        
        # Save as CSV
        df_matched.to_csv(args.output_path, index=False)
        
        print(f"\n✓ Completed in {time.time() - start_time:.1f}s")
        print(f"  Output: {args.output_path}")
        print(f"  Genes: {len(df_matched)}, Samples: {df.shape[1]}")
        print(f"  Unmatched regions: {unmatched_count}")
    else:
        print(f"\n✗ Error: No gene matches found. Check reference CSV format.")
        sys.exit(1)

if __name__ == "__main__":
    main()