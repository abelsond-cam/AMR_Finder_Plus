#!/usr/bin/env python3
"""
Debug script to process upstream raw protein sequence parquet files through AMR Finder Plus.

This script:
- Takes downstream files (with embeddings, without AMR columns)
- Extracts genome names from downstream files
- Searches upstream raw protein sequence files for matching genomes
- Processes matching upstream files through AMR Finder Plus with detailed debugging
"""

import pandas as pd
from pathlib import Path
from datetime import datetime
import sys

# Add the amrf directory to the path so we can import
sys.path.insert(0, str(Path(__file__).parent / "amrf"))

from run_amrf_genome_file import process_parquet_file


def main():
    """Main function to process upstream parquet files with debugging"""
    
    # Configuration
    downstream_base_dir = "/home/dca36/rds/hpc-work/data/BacFormer/processed/ast_esm_embeddings/campylobacter_jejuni"
    upstream_base_dir = "/home/dca36/rds/hpc-work/data/BacFormer/raw/ast/protein_sequences/campylobacter_jejuni"
    n_files = 1  # Start with 1 file as requested
    
    # AMR Finder columns to check for
    amr_columns = ['element_symbol', 'element_name', 'amrf_type', 'amrf_subtype']
    
    print("="*80)
    print("AMR Finder Plus Debug Processing")
    print("="*80)
    print(f"Downstream directory: {downstream_base_dir}")
    print(f"Upstream directory: {upstream_base_dir}")
    print(f"Processing first {n_files} file(s) without AMR annotations")
    print()
    
    # Step 1: Find downstream parquet files WITHOUT AMR columns
    downstream_files = list(Path(downstream_base_dir).rglob("*.parquet"))
    print(f"Total downstream parquet files found: {len(downstream_files)}")
    
    files_without_amr = []
    for file_path in downstream_files:
        try:
            df = pd.read_parquet(file_path, engine='pyarrow')
            columns = df.columns.tolist()
            has_all_amr_columns = all(col in columns for col in amr_columns)
            if not has_all_amr_columns:
                files_without_amr.append(str(file_path))
        except Exception as e:
            print(f"Warning: Could not check {file_path}: {e}")
    
    print(f"Downstream files WITHOUT AMR columns: {len(files_without_amr)}")
    print()
    
    # Take first N downstream files to process
    downstream_files_to_process = files_without_amr[:n_files]
    
    if not downstream_files_to_process:
        print("ERROR: No downstream files found to process!")
        return
    
    # Step 2: Find upstream parquet files
    upstream_files = list(Path(upstream_base_dir).rglob("*.parquet"))
    print(f"Upstream parquet files found: {len(upstream_files)}")
    
    if not upstream_files:
        print("ERROR: No upstream parquet files found!")
        return
    
    # Step 3: For each downstream file, find matching upstream file(s)
    matching_upstream_files = []
    
    for downstream_file_path in downstream_files_to_process:
        print("\n" + "="*80)
        print("ANALYZING DOWNSTREAM FILE")
        print("="*80)
        print(f"File: {Path(downstream_file_path).name}")
        print(f"Full path: {downstream_file_path}")
        
        # Extract genome names from downstream file
        try:
            downstream_df = pd.read_parquet(downstream_file_path, engine='pyarrow')
            downstream_genome_names = set(downstream_df['genome_name'].unique())
            print(f"Genomes in downstream file: {len(downstream_genome_names)}")
            
            # Show first few genome names
            print("Sample genome names:", list(downstream_genome_names)[:5])
            
        except Exception as e:
            print(f"ERROR reading downstream file: {e}")
            continue
        
        # Search upstream files for matching genomes
        print("\nSearching upstream files for matching genomes...")
        for upstream_file_path in upstream_files:
            try:
                upstream_df = pd.read_parquet(upstream_file_path, engine='pyarrow')
                upstream_genome_names = set(upstream_df['genome_name'].unique())
                
                # Check if there's any overlap
                overlap = downstream_genome_names & upstream_genome_names
                
                if overlap:
                    print(f"\n  MATCH FOUND: {Path(upstream_file_path).name}")
                    print(f"    Upstream genomes: {len(upstream_genome_names)}")
                    print(f"    Overlapping genomes: {len(overlap)}")
                    print(f"    Match percentage: {len(overlap)/len(downstream_genome_names)*100:.1f}%")
                    
                    # Add to matching files list
                    matching_upstream_files.append({
                        'upstream_file': str(upstream_file_path),
                        'downstream_file': Path(downstream_file_path).name,
                        'downstream_genomes': len(downstream_genome_names),
                        'upstream_genomes': len(upstream_genome_names),
                        'overlap': len(overlap)
                    })
                    
                    # Stop after first match (or could continue to find all matches)
                    break
                    
            except Exception as e:
                print(f"  Warning: Could not read {upstream_file_path.name}: {e}")
                continue
    
    if not matching_upstream_files:
        print("\nERROR: No matching upstream files found!")
        return
    
    print("\n" + "="*80)
    print("MATCHES FOUND")
    print("="*80)
    for match in matching_upstream_files:
        print(f"\nDownstream: {match['downstream_file']}")
        print(f"Upstream: {Path(match['upstream_file']).name}")
        print(f"  Downstream genomes: {match['downstream_genomes']}")
        print(f"  Upstream genomes: {match['upstream_genomes']}")
        print(f"  Overlap: {match['overlap']}")
    
    # Step 4: Process matching upstream files
    # Create log file with timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = Path(f"debug_amr_processing_{timestamp}.log")
    
    print(f"\nLog file: {log_file}")
    print()
    
    all_stats = []
    
    for match_idx, match in enumerate(matching_upstream_files, 1):
        upstream_file_path = match['upstream_file']
        
        print("\n" + "="*80)
        print(f"PROCESSING UPSTREAM FILE {match_idx}/{len(matching_upstream_files)}")
        print("="*80)
        print(f"File: {Path(upstream_file_path).name}")
        print(f"Full path: {upstream_file_path}")
        print(f"Source downstream file: {match['downstream_file']}")
        print()
        
        # Log to both console and file
        with open(log_file, 'a') as f:
            f.write(f"\n{'='*80}\n")
            f.write(f"PROCESSING UPSTREAM FILE {match_idx}/{len(matching_upstream_files)}\n")
            f.write(f"{'='*80}\n")
            f.write(f"File: {Path(upstream_file_path).name}\n")
            f.write(f"Full path: {upstream_file_path}\n")
            f.write(f"Source downstream file: {match['downstream_file']}\n")
            f.write(f"Started at: {datetime.now()}\n")
            f.write(f"{'='*80}\n\n")
        
        try:
            # Process the upstream file with debug_verbose=True
            output_dir, stats = process_parquet_file(
                input_file=upstream_file_path,
                verbose=True,
                log_file=log_file,
                debug_verbose=True
            )
            
            all_stats.append(stats)
            
            print("\n" + "="*80)
            print("FILE PROCESSING COMPLETE")
            print("="*80)
            print(f"File: {Path(upstream_file_path).name}")
            print(f"  Genomes processed: {stats['genomes_processed']}")
            print(f"  Genomes skipped: {stats['genomes_skipped']}")
            print(f"  Unique proteins: {stats['unique_proteins']}")
            print(f"  Duplicate proteins: {stats['duplicate_proteins']}")
            print(f"  AMR hits found: {stats['amr_hits']}")
            print(f"  AMR results saved to: {output_dir}")
            print()
            
        except Exception as e:
            error_msg = f"ERROR processing {upstream_file_path}: {e}"
            print(error_msg)
            with open(log_file, 'a') as f:
                f.write(f"\n{error_msg}\n")
                import traceback
                f.write(traceback.format_exc())
    
    # Overall summary
    print("\n" + "="*80)
    print("OVERALL SUMMARY")
    print("="*80)
    
    if all_stats:
        total_genomes = sum(s['genomes_processed'] for s in all_stats)
        total_hits = sum(s['amr_hits'] for s in all_stats)
        total_proteins = sum(s['unique_proteins'] for s in all_stats)
        
        print(f"Upstream files processed: {len(all_stats)}")
        print(f"Total genomes processed: {total_genomes}")
        print(f"Total unique proteins: {total_proteins}")
        print(f"Total AMR hits found: {total_hits}")
        print(f"Average AMR hits per genome: {total_hits/total_genomes:.2f}" if total_genomes > 0 else "N/A")
        
        # Write summary to log file
        with open(log_file, 'a') as f:
            f.write(f"\n{'='*80}\n")
            f.write("OVERALL SUMMARY\n")
            f.write(f"{'='*80}\n")
            f.write(f"Upstream files processed: {len(all_stats)}\n")
            f.write(f"Total genomes processed: {total_genomes}\n")
            f.write(f"Total unique proteins: {total_proteins}\n")
            f.write(f"Total AMR hits found: {total_hits}\n")
            f.write(f"Average AMR hits per genome: {total_hits/total_genomes:.2f}\n" if total_genomes > 0 else "N/A\n")
    
    print(f"\nDetailed log saved to: {log_file}")
    print("\n" + "="*80)
    print("PROCESSING COMPLETE")
    print("="*80)


if __name__ == "__main__":
    main()
