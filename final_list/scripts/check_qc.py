#!/usr/bin/env python3

def check_qc(flagstat_file, complexity_file, min_mapping_rate=70, min_complexity=0.7):
    """Check QC metrics from flagstat and complexity files."""
    try:
        # Check mapping rate
        with open(flagstat_file) as f:
            for line in f:
                if "mapped (" in line:
                    mapping_rate = float(line.split("(")[1].split("%")[0].strip())
                    if mapping_rate < min_mapping_rate:
                        print(f"Failed mapping rate check: {mapping_rate}% < {min_mapping_rate}%")
                        return False

        # Check library complexity
        with open(complexity_file) as f:
            next(f)  # Skip header
            complexity = float(next(f).split("\t")[1])
            if complexity < min_complexity:
                print(f"Failed complexity check: {complexity} < {min_complexity}")
                return False

        return True
    except Exception as e:
        print(f"Error in QC check: {str(e)}")
        return False

if __name__ == "__main__":
    import sys
    try:
        if check_qc(snakemake.input.flagstat, 
                    snakemake.input.complexity,
                    snakemake.params.min_mapping_rate, 
                    snakemake.params.min_complexity):
            with open(snakemake.output[0], 'w') as f:
                f.write("PASS\n")
        else:
            sys.exit(1)
    except Exception as e:
        print(f"Script failed: {str(e)}")
        sys.exit(1)