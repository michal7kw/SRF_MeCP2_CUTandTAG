import pandas as pd
import os
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Split MeCP2 enrichment analysis files')
    parser.add_argument('--work-dir', required=True,
                      help='Directory containing the input files')
    return parser.parse_args()

def main():
    args = parse_args()
    os.chdir(args.work_dir)

    # Read the CSV file
    df = pd.read_csv('cpg_enrichment_parallel.csv')

    # Create three dataframes based on binding_type
    binding_types_df = {
        'both': df[df['binding_type'] == 'both'],
        'nsc_only': df[df['binding_type'] == 'nsc_only'],
        'neu_only': df[df['binding_type'] == 'neu_only']
    }

    print("analyze_mecp2_cpg_enrichment:")
    print(f"Total number of peaks: {len(df)}")
    # Save each dataframe
    for binding_type, df_subset in binding_types_df.items():
        output_file = f"cpg_enrichment_{binding_type}.csv"
        df_subset.to_csv(output_file, index=False)
        print(f"Saved {binding_type} peaks ({len(df_subset)} rows) to {output_file}")

if __name__ == '__main__':
    main()