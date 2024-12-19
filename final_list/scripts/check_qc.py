#!/usr/bin/env python3

import sys
import re
import pandas as pd
import numpy as np
import os

def parse_flagstat(flagstat_file):
    """Parse samtools flagstat output to get mapping rate."""
    try:
        with open(flagstat_file) as f:
            content = f.read()
            
        # First try to find the percentage in the mapped line
        for line in content.split('\n'):
            if 'mapped (' in line:
                match = re.search(r'\(([\d.]+)%', line)
                if match:
                    return float(match.group(1))
        
        # If no percentage found, calculate from raw numbers
        lines = content.split('\n')
        total = int(lines[0].split()[0])  # First number in first line
        
        # Find mapped reads line
        for line in lines:
            if 'mapped (' in line:
                mapped = int(line.split()[0])
                return (mapped / total * 100) if total > 0 else 0.0
                
        raise ValueError("Could not find mapping statistics")
        
    except Exception as e:
        print(f"Error parsing flagstat file {flagstat_file}: {str(e)}", file=sys.stderr)
        print("File contents:", file=sys.stderr)
        with open(flagstat_file) as f:
            print(f.read(), file=sys.stderr)
        return 0.0

def parse_complexity(complexity_file):
    """Parse preseq output to estimate library complexity."""
    try:
        # First try reading with pandas
        try:
            df = pd.read_csv(complexity_file, sep='\t', comment='#')
            if len(df) == 0:
                raise ValueError("Empty dataframe")
        except Exception as e:
            print(f"Pandas reading failed: {str(e)}", file=sys.stderr)
            # Fall back to manual parsing
            with open(complexity_file) as f:
                lines = [line.strip() for line in f if not line.startswith('#')]
                if not lines:
                    raise ValueError("No data lines found")
                    
                # Parse header
                header = lines[0].split('\t')
                if len(header) < 2:
                    raise ValueError(f"Invalid header: {header}")
                    
                # Parse data lines
                data = []
                for line in lines[1:]:
                    try:
                        values = line.split('\t')
                        if len(values) >= 2:
                            total = float(values[0])
                            distinct = float(values[1])
                            data.append((total, distinct))
                    except ValueError:
                        continue
                
                if not data:
                    raise ValueError("No valid data points found")
                    
                # Convert to dataframe
                df = pd.DataFrame(data, columns=['TOTAL_READS', 'EXPECTED_DISTINCT'])
        
        # Find row closest to 10M reads
        target_reads = 10000000.0
        closest_idx = (df['TOTAL_READS'] - target_reads).abs().idxmin()
        row = df.iloc[closest_idx]
        
        total_reads = float(row['TOTAL_READS'])
        distinct_reads = float(row['EXPECTED_DISTINCT'])
        
        complexity = distinct_reads/total_reads if total_reads > 0 else 0.0
        print(f"Complexity calculation for {complexity_file}:", file=sys.stderr)
        print(f"Total reads: {total_reads}", file=sys.stderr)
        print(f"Distinct reads: {distinct_reads}", file=sys.stderr)
        print(f"Complexity: {complexity}", file=sys.stderr)
        
        return complexity
            
    except Exception as e:
        print(f"Error parsing complexity file {complexity_file}: {str(e)}", file=sys.stderr)
        print("File contents:", file=sys.stderr)
        with open(complexity_file) as f:
            print(f.read()[:1000], file=sys.stderr)
        return 0.0

def main(flagstat_file, complexity_file, qc_pass_file, min_mapping_rate=70, min_complexity=0.7):
    """Check QC metrics and create pass/fail file."""
    try:
        # Check if input files exist
        for f in [flagstat_file, complexity_file]:
            if not os.path.exists(f):
                raise FileNotFoundError(f"Input file not found: {f}")
        
        # Get QC metrics
        mapping_rate = parse_flagstat(flagstat_file)
        print(f"Mapping rate for {flagstat_file}: {mapping_rate}%", file=sys.stderr)
        
        complexity = parse_complexity(complexity_file)
        print(f"Library complexity for {complexity_file}: {complexity}", file=sys.stderr)
        
        # Check if passes QC
        passes_mapping = mapping_rate >= min_mapping_rate
        passes_complexity = complexity >= min_complexity
        passes_qc = passes_mapping and passes_complexity
        
        # Write results
        with open(qc_pass_file, 'w') as f:
            f.write(f"Sample QC Report\n")
            f.write(f"================\n")
            f.write(f"Mapping rate: {mapping_rate:.2f}% (threshold: {min_mapping_rate}%)\n")
            f.write(f"Library complexity: {complexity:.3f} (threshold: {min_complexity})\n")
            f.write(f"QC Status: {'PASS' if passes_qc else 'FAIL'}\n")
            if not passes_qc:
                f.write("\nFailed checks:\n")
                if not passes_mapping:
                    f.write(f"- Mapping rate too low: {mapping_rate:.2f}% < {min_mapping_rate}%\n")
                if not passes_complexity:
                    f.write(f"- Library complexity too low: {complexity:.3f} < {min_complexity}\n")
            
        # Don't exit with error, just report the status
        print(f"QC {'passed' if passes_qc else 'failed'} but continuing pipeline", file=sys.stderr)
        return 0
        
    except Exception as e:
        print(f"Error during QC check: {str(e)}", file=sys.stderr)
        with open(qc_pass_file, 'w') as f:
            f.write(f"Error during QC check: {str(e)}\n")
        return 1

if __name__ == "__main__":
    # Get file paths from snakemake
    flagstat_file = snakemake.input.flagstat
    complexity_file = snakemake.input.complexity
    qc_pass_file = snakemake.output.qc_pass
    
    # Get parameters
    min_mapping_rate = snakemake.params.min_mapping_rate
    min_complexity = snakemake.params.min_complexity
    
    # Run main and exit with its return code
    sys.exit(main(flagstat_file, complexity_file, qc_pass_file, 
                 min_mapping_rate=min_mapping_rate, 
                 min_complexity=min_complexity))