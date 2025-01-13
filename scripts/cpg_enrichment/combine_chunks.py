import pandas as pd
import glob
import os
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Combine chunk files from CPG enrichment analysis')
    parser.add_argument('--results-dir', required=True,
                      help='Directory containing the chunk files')
    return parser.parse_args()

def read_chunk_safely(file_path):
    try:
        if os.path.getsize(file_path) == 0:
            print(f"Warning: {file_path} is empty", file=sys.stderr)
            return None
            
        df = pd.read_csv(file_path)
        if df.empty:
            print(f"Warning: {file_path} contains no data", file=sys.stderr)
            return None
            
        print(f"Successfully read {file_path} with {len(df)} rows")
        return df
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}", file=sys.stderr)
        return None

def main():
    args = parse_args()
    chunks = sorted(glob.glob(os.path.join(args.results_dir, "chunk_*.csv")))
    print(f"Found {len(chunks)} chunk files")

    dfs = []
    for chunk_file in chunks:
        df = read_chunk_safely(chunk_file)
        if df is not None:
            dfs.append(df)

    if not dfs:
        print("Error: No valid data found in any chunk files", file=sys.stderr)
        sys.exit(1)

    combined = pd.concat(dfs, ignore_index=True)
    output_file = os.path.join(args.results_dir, "mecp2_cpg_enrichment_parallel.csv")
    combined.to_csv(output_file, index=False)
    print(f"Combined {len(dfs)} valid chunks into {output_file}")
    print(f"Total rows in combined file: {len(combined)}")

    for chunk_file, df in zip(chunks, dfs):
        print(f"Chunk {os.path.basename(chunk_file)}: {len(df)} rows")

if __name__ == '__main__':
    main() 