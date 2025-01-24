# Cross Analysis Pipeline

This pipeline analyzes methylation patterns and SMARCB1 binding in different gene sets.

## Dependencies

Required Python packages:
- pandas
- numpy
- pyBigWig
- matplotlib
- seaborn
- scipy
- pybedtools

## Directory Structure

```
Cross_final/
├── src/
│   ├── config.py
│   ├── methylation_analysis.py
│   └── genomic_utils.py
├── data/
│   ├── up.csv
│   ├── down.csv
│   └── no_deg.csv
└── results/
```

## Usage

1. Copy your gene list files (up.csv, down.csv, no_deg.csv) to the `data/` directory
2. Update paths in `config.py` if necessary
3. Run the analysis:
   ```bash
   python src/methylation_analysis.py
   ```

## Output

The pipeline generates:
1. Methylation profile plots for each condition
2. SMARCB1 binding enrichment analysis
3. Statistical test results for differences between gene categories
