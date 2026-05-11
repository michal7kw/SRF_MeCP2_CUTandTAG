import csv
from operator import itemgetter
from typing import List, Dict, Tuple
import argparse

def load_and_filter_cpg(csv_path: str, signal_threshold: float = 0.1) -> List[Dict]:
    """Load and filter CPG enrichment data"""
    filtered_data = []
    
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Convert string values to float
            row['exo_signal'] = float(row['exo_signal'])
            row['endo_signal'] = float(row['endo_signal'])
            row['enrichment'] = float(row['enrichment'])
            row['start'] = int(row['start'])
            row['end'] = int(row['end'])
            
            # Filter based on threshold
            if row['exo_signal'] > signal_threshold or row['endo_signal'] > signal_threshold:
                filtered_data.append(row)
    
    # Sort by enrichment score in descending order
    return sorted(filtered_data, key=lambda x: x['enrichment'], reverse=True)

def parse_gtf_attributes(attribute_string: str) -> Dict[str, str]:
    """Parse GTF attribute string into a dictionary"""
    attributes = {}
    for attribute in attribute_string.strip().split(';'):
        if attribute:
            try:
                key, value = attribute.strip().split(' ', 1)
                attributes[key] = value.strip('"')
            except ValueError:
                continue
    return attributes

def load_gtf(gtf_path: str) -> List[Dict]:
    """Load and process GTF file"""
    genes = []
    
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
                
            attributes = parse_gtf_attributes(fields[8])
            if 'gene_name' not in attributes or 'gene_type' not in attributes:
                continue
                
            genes.append({
                'seqname': fields[0],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'gene_name': attributes['gene_name'],
                'gene_type': attributes['gene_type']
            })
    
    return genes

def determine_genomic_region(peak: Dict, genes: List[Dict]) -> str:
    """Determine the genomic region for a given CPG site"""
    chr_genes = [g for g in genes if g['seqname'] == peak['chr']]
    
    # Check if the peak overlaps with any gene
    for gene in chr_genes:
        if (peak['start'] <= gene['end'] and peak['end'] >= gene['start']):
            return f"{gene['gene_name']} ({gene['gene_type']})"
    
    # If no overlap found, find nearest gene
    distances = []
    for gene in chr_genes:
        if peak['start'] > gene['end']:
            distance = peak['start'] - gene['end']
        elif peak['end'] < gene['start']:
            distance = gene['start'] - peak['end']
        else:
            distance = 0
        distances.append((distance, gene['gene_name'], gene['gene_type']))
    
    if distances:
        min_distance, nearest_gene, gene_type = min(distances, key=itemgetter(0))
        if min_distance > 10000:  # More than 10kb away
            return "Intergenic"
        else:
            return f"Near {nearest_gene} ({gene_type}), {min_distance}bp"
    
    return "Intergenic"

def write_results(data: List[Dict], output_file: str):
    """Write results to CSV file"""
    if not data:
        return
        
    fieldnames = list(data[0].keys())
    
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Process CpG enrichment data')
    parser.add_argument('--working-dir', type=str, required=True,
                       help='Path to working directory')
    parser.add_argument('--data-dir', type=str, required=True,
                       help='Path to data directory')
    parser.add_argument('--results-dir', type=str, required=True,
                       help='Path to results directory')
    args = parser.parse_args()

    # File paths
    cpg_file = f"{args.results_dir}/cpg_enrichment_NSC.csv"
    gtf_file = f"{args.data_dir}/gencode.vM10.annotation.gtf"
    output_file = f"{args.results_dir}/cpg_enrichment_annotated.csv"
    
    # Load and filter CPG data
    print("Loading and filtering CPG data...")
    cpg_data = load_and_filter_cpg(cpg_file)
    
    # Load GTF data
    print("Loading GTF data...")
    genes = load_gtf(gtf_file)
    
    # Add genomic region annotation
    print("Annotating genomic regions...")
    for peak in cpg_data:
        peak['genomic_region'] = determine_genomic_region(peak, genes)
    
    # Save results
    write_results(cpg_data, output_file)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main() 