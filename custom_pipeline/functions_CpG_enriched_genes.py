# Standard library imports
import os
import subprocess
from collections import defaultdict
from pathlib import Path

# Third party imports
import numpy as np
import pandas as pd
from scipy import stats
import statistics

# Visualization
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn3
import glob
from IPython.display import Image, display
from venn import venn

def get_peaks_with_cpg(peak_file, cpg_file, extend=300, coverage_threshold=20):
    """
    Identifies peaks that overlap with CpG islands above a coverage threshold.
    
    Args:
        peak_file: BED file containing peak regions
        cpg_file: BED file containing CpG islands
        extend: Number of base pairs to extend CpG islands (default: 300bp)
        coverage_threshold: Minimum percentage overlap required (default: 20%)
    
    Returns:
        dict: Peak coordinates that meet the CpG overlap criteria
    """
    print(f"\nProcessing {peak_file}")
    
    # Validate input files
    if not all(os.path.exists(f) for f in [peak_file, cpg_file]):
        print("Error: Input file(s) missing!")
        return set()
    
    # Filter for standard chromosomes (chr1-22, X, Y)
    standard_chroms = [f"chr{i}" for i in range(1, 23)] + ['chrX', 'chrY']
    
    # Create temporary filtered files
    filtered_peaks = "temp_filtered_peaks.bed"
    filtered_cpg = "temp_filtered_cpg.bed"
    
    # Filter peaks and count
    with open(peak_file) as fin, open(filtered_peaks, 'w') as fout:
        peaks = [line for line in fin if line.split('\t')[0] in standard_chroms]
        fout.writelines(peaks)
        print(f"Peaks after chromosome filtering: {len(peaks)}")
    
    # Filter and extend CpG islands
    with open(cpg_file) as fin, open(filtered_cpg, 'w') as fout:
        fout.writelines(line for line in fin if line.split('\t')[0] in standard_chroms)
    
    # # Print first 5 lines from each file for debugging
    # with open(filtered_peaks) as f:
    #     print("\nFirst 5 lines of filtered peaks:")
    #     for i, line in enumerate(f):
    #         if i < 5:
    #             print(line.strip())  
    # with open(filtered_cpg) as f:
    #     print("\nFirst 5 lines of filtered CpG:")
    #     for i, line in enumerate(f):
    #         if i < 5:
    #             print(line.strip())
    
    # Extend CpG regions and find overlaps using bedtools
    subprocess.run(f"bedtools slop -i {filtered_cpg} -g DATA/genome.size -b {extend} > temp_extended_cpg.bed", shell=True)
    result = subprocess.run(f"bedtools intersect -a {filtered_peaks} -b temp_extended_cpg.bed -wao", 
                          shell=True, capture_output=True, text=True)
    
    # # Debug: Print first few lines of bedtools output
    # print("\nFirst 5 lines of bedtools intersect output:")
    # for i, line in enumerate(result.stdout.strip().split('\n')[:5]):
    #     print(line)
    
    # Process overlaps and calculate coverage
    peaks_with_cpg = {}
    coverage_stats = []
    current_peak = {'id': None, 'overlaps': []}
    qualified_peaks = 0
    total_peaks_with_overlap = 0
    
    for line in result.stdout.strip().split('\n'):
        fields = line.split('\t')
        
        # Fields from peak (first BED file)
        peak_chr = fields[0]
        peak_start = fields[1]
        peak_end = fields[2]
        peak_id = f"{peak_chr}:{peak_start}-{peak_end}"
        
        # Check if there's a valid CpG overlap
        if len(fields) > 6 and fields[4] != "." and fields[4] != "-1":  # Changed condition
            try:
                cpg_start = int(fields[5])  # Changed from fields[6]
                cpg_end = int(fields[6])    # Changed from fields[7]
                overlap_amount = int(fields[-1])  # Last field is overlap amount
                
                # If we encounter a new peak
                if peak_id != current_peak['id']:
                    # Process previous peak if it exists
                    if current_peak['overlaps']:
                        max_coverage = max(current_peak['overlaps'])
                        total_peaks_with_overlap += 1
                        if max_coverage >= coverage_threshold:
                            coverage_stats.append(max_coverage)
                            peaks_with_cpg[current_peak['id']] = max_coverage
                            qualified_peaks += 1
                    
                    current_peak = {'id': peak_id, 'overlaps': []}
                
                # Calculate overlap coverage using CpG length
                if overlap_amount > 0:  # If there is overlap
                    cpg_length = cpg_end - cpg_start  # CpG island length
                    if cpg_length > 0:  # Avoid division by zero
                        coverage = (overlap_amount / cpg_length) * 100
                        current_peak['overlaps'].append(coverage)
                        
                        # # Debug: Print coverage calculation for first few peaks
                        # if qualified_peaks < 3:
                        #     print(f"\nCoverage calculation for {peak_id}:")
                        #     print(f"CpG region: {cpg_start}-{cpg_end}")
                        #     print(f"Overlap amount: {overlap_amount}")
                        #     print(f"CpG length: {cpg_length}")
                        #     print(f"Coverage: {coverage:.2f}%")
            
            except (ValueError, IndexError) as e:
                continue
        
        else:
            # If we encounter a new peak with no overlap
            if peak_id != current_peak['id']:
                # Process previous peak if it exists
                if current_peak['overlaps']:
                    max_coverage = max(current_peak['overlaps'])
                    total_peaks_with_overlap += 1
                    if max_coverage >= coverage_threshold:
                        coverage_stats.append(max_coverage)
                        peaks_with_cpg[current_peak['id']] = max_coverage
                        qualified_peaks += 1
                
                current_peak = {'id': peak_id, 'overlaps': []}
    
    # Process the last peak
    if current_peak['overlaps']:
        max_coverage = max(current_peak['overlaps'])
        total_peaks_with_overlap += 1
        if max_coverage >= coverage_threshold:
            coverage_stats.append(max_coverage)
            peaks_with_cpg[current_peak['id']] = max_coverage
            qualified_peaks += 1
    
    # Print coverage statistics
    if coverage_stats:
        print(f"\nPeaks with any CpG overlap: {total_peaks_with_overlap}")
        print(f"Peaks meeting {coverage_threshold}% coverage threshold: {qualified_peaks}")
        print(f"Coverage range: {min(coverage_stats):.2f}% - {max(coverage_stats):.2f}%")
        print(f"Mean coverage of qualified peaks: {statistics.mean(coverage_stats):.2f}%")
    else:
        print("\nNo peaks met the coverage criteria")
    
    # Cleanup temporary files
    for f in [filtered_peaks, filtered_cpg, "temp_extended_cpg.bed"]:
        if os.path.exists(f):
            os.remove(f)
    
    return peaks_with_cpg

def get_genes_with_cpg_enrichment(peak_file, cpg_file, gtf_file, output_dir, cell_type, condition, 
                                 extend_cpg=300, extend_tss=2000, coverage_threshold=20):
    """
    Identifies genes with CpG-overlapping peaks near their TSS regions.
    """
    print(f"\nAnalyzing {cell_type} {condition} peaks")
    
    # Validate input files
    if not os.path.exists(gtf_file):
        print(f"Error: GTF file not found: {gtf_file}")
        return None
        
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Get peaks overlapping with CpG islands
    cpg_peaks = get_peaks_with_cpg(peak_file, cpg_file, extend=extend_cpg, 
                                  coverage_threshold=coverage_threshold)
    
    if not cpg_peaks:
        print(f"No CpG-overlapping peaks found for {cell_type} {condition}")
        return {'genes': pd.DataFrame(), 'total_peaks': 0}
    
    # 2. Extract and extend TSS regions
    try:
        tss_bed = os.path.join(output_dir, "temp_tss.bed")
        extract_tss_regions(gtf_file, tss_bed)
        
        extended_tss = os.path.join(output_dir, "temp_extended_tss.bed")
        subprocess.run(f"bedtools slop -i {tss_bed} -g DATA/genome.size -b {extend_tss} > {extended_tss}",
                      shell=True, check=True)
        
        # 3. Write CpG-overlapping peaks to temporary file
        temp_peaks = os.path.join(output_dir, "temp_cpg_peaks.bed")
        with open(temp_peaks, 'w') as f:
            for peak, coverage in cpg_peaks.items():
                chrom, pos = peak.split(':')
                start, end = pos.split('-')
                f.write(f"{chrom}\t{start}\t{end}\t{coverage}\n")
        
        # 4. Find overlaps between peaks and TSS regions
        result = subprocess.run(
            f"bedtools intersect -a {extended_tss} -b {temp_peaks} -wo",
            shell=True, capture_output=True, text=True, check=True
        )
        
        # 5. Process results
        gene_data = defaultdict(lambda: {
            'peaks': set(),
            'total_peak_coverage': 0,
            'distance_to_tss': [],
            'peak_coverages': []
        })
        
        for line in result.stdout.strip().split('\n'):
            if not line:
                continue
            
            fields = line.split('\t')
            # print(f"fields: {fields}")
            gene_name = fields[3]
            tss_pos = int(fields[1]) + extend_tss
            peak_start = int(fields[5])
            peak_end = int(fields[6])
            peak_id = f"{fields[4]}:{fields[5]}-{fields[6]}"
            peak_coverage = float(fields[7])
            
            peak_center = (peak_start + peak_end) // 2
            distance = abs(peak_center - tss_pos)
            
            gene_data[gene_name]['peaks'].add(peak_id)
            gene_data[gene_name]['total_peak_coverage'] += int(fields[-1])
            gene_data[gene_name]['distance_to_tss'].append(distance)
            gene_data[gene_name]['peak_coverages'].append(peak_coverage)
        
        # 6. Create summary DataFrame
        summary_data = []
        for gene, data in gene_data.items():
            if data['peaks']:  # Only include genes with peaks
                summary_data.append({
                    'gene_name': gene,
                    'num_peaks': len(data['peaks']),
                    'total_coverage': data['total_peak_coverage'],
                    'avg_peak_size': data['total_peak_coverage'] / len(data['peaks']),
                    'min_distance_to_tss': min(data['distance_to_tss']),
                    'mean_cpg_coverage': statistics.mean(data['peak_coverages']),
                    'peaks': ','.join(data['peaks'])
                })
        
        df = pd.DataFrame(summary_data)
        
        # 7. Save results
        if not df.empty:
            output_file = os.path.join(output_dir, f"{cell_type}_{condition}_cpg_genes.tsv")
            df.sort_values('num_peaks', ascending=False).to_csv(output_file, sep='\t', index=False)
            print(f"Found {len(df)} genes with CpG-overlapping peaks")
            print(f"Results saved to: {output_file}")
        else:
            print(f"No genes found with CpG-overlapping peaks for {cell_type} {condition}")
        
        return {'genes': df, 'total_peaks': len(cpg_peaks)}
        
    except Exception as e:
        print(f"Error processing gene enrichment: {str(e)}")
        return {'genes': pd.DataFrame(), 'total_peaks': 0}
    finally:
        # Cleanup temporary files
        for f in [tss_bed, extended_tss, temp_peaks]:
            if os.path.exists(f):
                os.remove(f)

def extract_tss_regions(gtf_file, output_bed):
    """Extracts TSS regions from GTF file."""
    with open(gtf_file) as fin, open(output_bed, 'w') as fout:
        for line in fin:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if fields[2] != "gene":
                continue
            
            # Parse gene information
            info_dict = {}
            for item in fields[8].rstrip(';').split('; '):
                try:
                    key, value = item.strip().split(' ', 1)
                    info_dict[key] = value.strip('"')
                except ValueError:
                    continue
            
            gene_name = info_dict.get('gene_name', '')
            gene_type = info_dict.get('gene_type', '')
            
            if not gene_name or gene_type != "protein_coding":
                continue
            
            # Get TSS position
            chrom = fields[0]
            if fields[6] == '+':
                tss_pos = int(fields[3])
                fout.write(f"{chrom}\t{tss_pos-1}\t{tss_pos}\t{gene_name}\n")
            else:
                tss_pos = int(fields[4])
                fout.write(f"{chrom}\t{tss_pos-1}\t{tss_pos}\t{gene_name}\n")

def create_comparison_summary(results, output_dir):
    """Creates a summary of gene overlaps between conditions."""
    summary_file = os.path.join(output_dir, "comparison_summary.txt")
    
    with open(summary_file, 'w') as f:
        for cell_type in results:
            exo_genes = set(results[cell_type]['Exo']['genes']['gene_name']) if not results[cell_type]['Exo']['genes'].empty else set()
            endo_genes = set(results[cell_type]['Endo']['genes']['gene_name']) if not results[cell_type]['Endo']['genes'].empty else set()
            
            common_genes = exo_genes & endo_genes
            
            f.write(f"\n{cell_type} Analysis:\n")
            f.write(f"Exogenous-specific genes: {len(exo_genes - common_genes)}\n")
            f.write(f"Endogenous-specific genes: {len(endo_genes - common_genes)}\n")
            f.write(f"Common genes: {len(common_genes)}\n")
            
            if common_genes:
                common_genes_file = os.path.join(output_dir, f"{cell_type}_common_genes.txt")
                with open(common_genes_file, 'w') as cf:
                    for gene in sorted(common_genes):
                        cf.write(f"{gene}\n")

def load_and_process_data(results_dir):
    """
    Load and process data with CpG coverage calculations.
    """
    cell_types = ['NSC', 'Neuron']
    conditions = ['Exo', 'Endo']
    data = {}
    
    for cell_type in cell_types:
        data[cell_type] = {}
        for condition in conditions:
            file_path = os.path.join(results_dir, f"{cell_type}_{condition}_cpg_genes.tsv")
            if os.path.exists(file_path):
                df = pd.read_csv(file_path, sep='\t')
                data[cell_type][condition] = df
            else:
                print(f"Warning: Missing file {file_path}")
                data[cell_type][condition] = pd.DataFrame()
    
    return data

def compare_endo_cpg_coverage(data):
    """
    Compare CpG islands coverage of Endo MeCP2 between NPCs and Neurons.
    """
    # Extract Endo data for both cell types
    nsc_endo = data['NSC']['Endo']
    neuron_endo = data['Neuron']['Endo']
    
    # Create comparison plots
    plt.figure(figsize=(12, 6))
    
    # Plot 1: Coverage Distribution
    plt.subplot(1, 2, 1)
    sns.kdeplot(data=nsc_endo['mean_cpg_coverage'], label='NSCs', color='#bcdae6')
    sns.kdeplot(data=neuron_endo['mean_cpg_coverage'], label='Neurons', color='#efb3b1')
    plt.xlabel('Mean CpG Coverage (%)')
    plt.ylabel('Density')
    plt.title('Endo MeCP2 CpG Coverage Distribution')
    plt.legend()
    
    # Plot 2: Box Plot
    plt.subplot(1, 2, 2)
    plot_data = pd.concat([
        pd.DataFrame({'Cell Type': 'NSCs', 'Coverage': nsc_endo['mean_cpg_coverage']}),
        pd.DataFrame({'Cell Type': 'Neurons', 'Coverage': neuron_endo['mean_cpg_coverage']})
    ])
    sns.boxplot(data=plot_data, x='Cell Type', y='Coverage')
    plt.title('CpG Coverage Comparison')
    
    plt.tight_layout()
    plt.show()
    
    # Statistical comparison
    stat_result = stats.mannwhitneyu(nsc_endo['mean_cpg_coverage'],
                                   neuron_endo['mean_cpg_coverage'],
                                   alternative='two-sided')
    
    # Display statistical results
    print("Mann-Whitney U test results:")
    print(f"Statistic: {stat_result.statistic}")
    print(f"p-value: {stat_result.pvalue}")
    print("\nSummary Statistics:")
    print(f"NSCs mean coverage: {nsc_endo['mean_cpg_coverage'].mean():.2f}%")
    print(f"Neurons mean coverage: {neuron_endo['mean_cpg_coverage'].mean():.2f}%")
    print(f"NSCs median coverage: {nsc_endo['mean_cpg_coverage'].median():.2f}%")
    print(f"Neurons median coverage: {neuron_endo['mean_cpg_coverage'].median():.2f}%")

    # Create identical plots for Exo comparison
    nsc_exo = data['NSC']['Exo']
    neuron_exo = data['Neuron']['Exo']
    
    plt.figure(figsize=(12, 6))
    
    # Plot 1: Coverage Distribution for Exo
    plt.subplot(1, 2, 1)
    sns.kdeplot(data=nsc_exo['mean_cpg_coverage'], label='NSCs', color='#1919f5')
    sns.kdeplot(data=neuron_exo['mean_cpg_coverage'], label='Neurons', color='#ec4b3d')
    plt.xlabel('Mean CpG Coverage (%)')
    plt.ylabel('Density')
    plt.title('Exo MeCP2 CpG Coverage Distribution')
    plt.legend()
    
    # Plot 2: Box Plot for Exo
    plt.subplot(1, 2, 2)
    plot_data = pd.concat([
        pd.DataFrame({'Cell Type': 'NSCs', 'Coverage': nsc_exo['mean_cpg_coverage']}),
        pd.DataFrame({'Cell Type': 'Neurons', 'Coverage': neuron_exo['mean_cpg_coverage']})
    ])
    sns.boxplot(data=plot_data, x='Cell Type', y='Coverage')
    plt.title('CpG Coverage Comparison')
    
    plt.tight_layout()
    plt.show()
    
    # Statistical comparison for Exo
    stat_result_exo = stats.mannwhitneyu(nsc_exo['mean_cpg_coverage'],
                                       neuron_exo['mean_cpg_coverage'],
                                       alternative='two-sided')
    
    # Display statistical results for Exo
    print("\nExo MeCP2 Mann-Whitney U test results:")
    print(f"Statistic: {stat_result_exo.statistic}")
    print(f"p-value: {stat_result_exo.pvalue}")
    print("\nSummary Statistics:")
    print(f"NSCs mean coverage: {nsc_exo['mean_cpg_coverage'].mean():.2f}%")
    print(f"Neurons mean coverage: {neuron_exo['mean_cpg_coverage'].mean():.2f}%")
    print(f"NSCs median coverage: {nsc_exo['mean_cpg_coverage'].median():.2f}%")
    print(f"Neurons median coverage: {neuron_exo['mean_cpg_coverage'].median():.2f}%")

def analyze_common_cpg_targets(data):
    """
    Analyze CpG targets common between all four experimental groups.
    """
    # Extract gene sets
    nsc_endo = set(data['NSC']['Endo']['gene_name'])
    neuron_endo = set(data['Neuron']['Endo']['gene_name'])
    nsc_exo = set(data['NSC']['Exo']['gene_name'])
    neuron_exo = set(data['Neuron']['Exo']['gene_name'])
    
    # Calculate common targets across all groups
    all_common = nsc_endo.intersection(neuron_endo, nsc_exo, neuron_exo)
    all_unique = nsc_endo.union(neuron_endo, nsc_exo, neuron_exo)
    
    # Calculate percentage of common targets
    common_percent = (len(all_common) / len(all_unique)) * 100
    
    # Create the 4-way Venn diagram
    plt.figure(figsize=(10, 10))
    venn({
        'NSC Endo': nsc_endo,
        'Neuron Endo': neuron_endo,
        'NSC Exo': nsc_exo,
        'Neuron Exo': neuron_exo
    })
    plt.title('CpG Targets Across All Groups')
    plt.tight_layout()
    plt.show()
    
    # Display overlap statistics
    print("CpG Target Statistics:")
    print("-" * 50)
    print(f"NSC Endo targets: {len(nsc_endo)}")
    print(f"Neuron Endo targets: {len(neuron_endo)}")
    print(f"NSC Exo targets: {len(nsc_exo)}")
    print(f"Neuron Exo targets: {len(neuron_exo)}")
    print("-" * 50)
    print(f"Total unique targets: {len(all_unique)}")
    print(f"Common targets (all groups): {len(all_common)} ({common_percent:.2f}%)")
    
    # Calculate pairwise overlaps
    pairs = [
        ('NSC Endo', 'Neuron Endo', nsc_endo, neuron_endo),
        ('NSC Endo', 'NSC Exo', nsc_endo, nsc_exo),
        ('NSC Endo', 'Neuron Exo', nsc_endo, neuron_exo),
        ('Neuron Endo', 'NSC Exo', neuron_endo, nsc_exo),
        ('Neuron Endo', 'Neuron Exo', neuron_endo, neuron_exo),
        ('NSC Exo', 'Neuron Exo', nsc_exo, neuron_exo)
    ]
    
    print("\nPairwise Overlaps:")
    print("-" * 50)
    for label1, label2, set1, set2 in pairs:
        overlap = len(set1.intersection(set2))
        total = len(set1.union(set2))
        percent = (overlap / total) * 100
        print(f"{label1} vs {label2}: {overlap} ({percent:.2f}%)")
    
    # Calculate group-specific targets
    print("\nGroup-Specific Targets:")
    print("-" * 50)
    nsc_endo_specific = len(nsc_endo - (neuron_endo | nsc_exo | neuron_exo))
    neuron_endo_specific = len(neuron_endo - (nsc_endo | nsc_exo | neuron_exo))
    nsc_exo_specific = len(nsc_exo - (nsc_endo | neuron_endo | neuron_exo))
    neuron_exo_specific = len(neuron_exo - (nsc_endo | neuron_endo | nsc_exo))
    
    print(f"NSC Endo-specific: {nsc_endo_specific}")
    print(f"Neuron Endo-specific: {neuron_endo_specific}")
    print(f"NSC Exo-specific: {nsc_exo_specific}")
    print(f"Neuron Exo-specific: {neuron_exo_specific}")

def analyze_exo_vs_endo_enrichment(data):
    """
    Analyze enrichment of Exo vs Endo MeCP2 CpG coverage.
    """
    enriched_genes = {'NSC': [], 'Neuron': []}
    
    for cell_type in ['NSC', 'Neuron']:
        exo_df = data[cell_type]['Exo']
        endo_df = data[cell_type]['Endo']
        
        # Merge Exo and Endo data
        merged = pd.merge(exo_df, endo_df, 
                         on='gene_name', 
                         suffixes=('_exo', '_endo'))
        
        # Calculate enrichment ratio using mean_cpg_coverage
        merged['enrichment_ratio'] = merged['mean_cpg_coverage_exo'] / merged['mean_cpg_coverage_endo']
        
        # Filter for significantly enriched genes (log2 fold change >= 0.5)
        merged['log2_enrichment'] = np.log2(merged['enrichment_ratio'])
        enriched = merged[merged['log2_enrichment'] >= 0.5].sort_values('enrichment_ratio', ascending=False)
        enriched_genes[cell_type] = enriched
        
        # Save enriched genes to TSV files
        output_file = f"results/enriched_genes_{cell_type}.tsv"
        enriched.to_csv(output_file, sep='\t', index=False)
        
        # Optional: Print summary statistics
        print(f"\n{cell_type} Enrichment Summary:")
        print(f"Total genes analyzed: {len(merged)}")
        print(f"Enriched genes: {len(enriched)}")
        print(f"Mean enrichment ratio: {merged['enrichment_ratio'].mean():.2f}")
        print(f"Median enrichment ratio: {merged['enrichment_ratio'].median():.2f}")
    
    # Find common enriched genes
    common_enriched = set(enriched_genes['NSC']['gene_name']).intersection(
        set(enriched_genes['Neuron']['gene_name']))
    
    print("\nGenes enriched in Exo vs Endo in both NPCs and Neurons:")
    for gene in sorted(common_enriched):
        print(gene)

def generate_comprehensive_analysis(results_dir):
    """
    Generate all analyses and visualizations.
    """
    # Load data
    data = load_and_process_data(results_dir)
    
    # Generate all analyses
    # compare_endo_cpg_coverage(data)
    analyze_common_cpg_targets(data)
    analyze_exo_vs_endo_enrichment(data)

def load_data(results_dir):
    """
    Load all result files into a dictionary of DataFrames.
    
    Args:
        results_dir (str): Path to the results directory
    
    Returns:
        dict: Nested dictionary containing DataFrames for each cell type and condition
    """
    cell_types = ['NSC', 'Neuron']
    conditions = ['Exo', 'Endo']
    data = {}
    
    for cell_type in cell_types:
        data[cell_type] = {}
        for condition in conditions:
            file_path = os.path.join(results_dir, f"{cell_type}_{condition}_cpg_genes.tsv")
            if os.path.exists(file_path):
                data[cell_type][condition] = pd.read_csv(file_path, sep='\t')
            else:
                print(f"Warning: Missing file {file_path}")
                data[cell_type][condition] = pd.DataFrame()
                
    return data

def plot_peak_distribution(data):
    """
    Create violin plots showing the distribution of peaks per gene.
    """
    plt.figure(figsize=(10, 6))
    data_points = []
    
    for cell_type in data:
        for condition in data[cell_type]:
            if not data[cell_type][condition].empty:
                peak_counts = data[cell_type][condition]['num_peaks']
                label = f"{cell_type}-{condition}"
                data_points.extend([(x, label) for x in peak_counts])
    
    if data_points:
        df = pd.DataFrame(data_points, columns=['Peaks', 'Group'])
        sns.violinplot(data=df, x='Group', y='Peaks')
        plt.xticks(rotation=45)
        plt.title('Distribution of CpG-overlapping Peaks per Gene')
        plt.tight_layout()
        plt.show()

def plot_distance_distribution(data):
    """
    Create KDE plots showing the distribution of distances to TSS.
    """
    plt.figure(figsize=(12, 6))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    color_idx = 0
    
    for cell_type in data:
        for condition in data[cell_type]:
            if not data[cell_type][condition].empty:
                distances = data[cell_type][condition]['min_distance_to_tss']
                sns.kdeplot(data=distances, label=f"{cell_type}-{condition}",
                          color=colors[color_idx])
                color_idx += 1
    
    plt.xlabel('Distance to TSS (bp)')
    plt.ylabel('Density')
    plt.title('Distribution of Minimum Distances to TSS')
    plt.legend()
    plt.tight_layout()
    plt.show()

def create_venn_diagrams(data):
    """
    Create Venn diagrams showing overlap between conditions for each cell type.
    """
    colors = {'NSC': ['#1f77b4', '#2ca02c'], 'Neuron': ['#ff7f0e', '#d62728']}
    
    for cell_type in data:
        plt.figure(figsize=(8, 8))
        
        exo_genes = set(data[cell_type]['Exo']['gene_name']) if not data[cell_type]['Exo'].empty else set()
        endo_genes = set(data[cell_type]['Endo']['gene_name']) if not data[cell_type]['Endo'].empty else set()
        
        # Create Venn diagram
        venn2([exo_genes, endo_genes],
              set_labels=('Exogenous', 'Endogenous'),
              set_colors=colors[cell_type])
        
        # Add title separately
        plt.title(f'{cell_type} Gene Overlap')
        
        plt.tight_layout()
        plt.show()

def create_venn_diagrams_percentages(data):
    """
    Create Venn diagrams showing overlap between conditions for each cell type.
    Percentages are calculated relative to the Exogenous sample size.
    """
    colors = {'NSC': ['#1919f5', '#bcdae6'], 'Neuron': ['#ff0000', '#efb3b1']}
    
    for cell_type in data:
        plt.figure(figsize=(8, 8))
        
        exo_genes = set(data[cell_type]['Exo']['gene_name']) if not data[cell_type]['Exo'].empty else set()
        endo_genes = set(data[cell_type]['Endo']['gene_name']) if not data[cell_type]['Endo'].empty else set()
        
        # Calculate total number of exo genes
        total_exo_genes = len(exo_genes)
        
        # Create Venn diagram with percentages relative to exo sample
        venn2([exo_genes, endo_genes],
              set_labels=('', ''),
              set_colors=colors[cell_type], alpha=0.5,
              subset_label_formatter=lambda x: f'{(x/total_exo_genes)*100:.1f}%' if total_exo_genes > 0 else '0%')
        
        # Add title separately
        plt.title(f'{cell_type} Gene Overlap (% of Exogenous Genes)')
        
        plt.tight_layout()
        plt.show()


def create_venn_diagrams_endo_comparison(data):
    """
    Create Venn diagram showing overlap between NSCs and Neurons for Endogenous MeCP2.
    """
    plt.figure(figsize=(8, 8))
    
    # Get Endogenous genes for each cell type
    nsc_endo = set(data['NSC']['Endo']['gene_name']) if not data['NSC']['Endo'].empty else set()
    neuron_endo = set(data['Neuron']['Endo']['gene_name']) if not data['Neuron']['Endo'].empty else set()
    
    # Calculate total number of genes across both sets
    total_genes = len(nsc_endo.union(neuron_endo))
    
    # Create Venn diagram with percentages relative to total genes
    venn2([nsc_endo, neuron_endo],
          set_labels=('', ''),
          set_colors=['#bcdae6', '#efb3b1'], alpha=0.8,
          subset_label_formatter=lambda x: f'{(x/total_genes)*100:.1f}%' if total_genes > 0 else '0%')
    
    # Add title
    plt.title('Endogenous MeCP2 Gene Overlap\n(% of Total Genes)')
    
    plt.tight_layout()
    plt.show()
    
    # Print statistics
    common_genes = nsc_endo.intersection(neuron_endo)
    nsc_specific = nsc_endo - neuron_endo
    neuron_specific = neuron_endo - nsc_endo
    
    print("\nEndogenous MeCP2 Binding Statistics:")
    print(f"NSC-specific genes: {len(nsc_specific)} ({len(nsc_specific)/total_genes*100:.1f}% of total genes)")
    print(f"Neuron-specific genes: {len(neuron_specific)} ({len(neuron_specific)/total_genes*100:.1f}% of total genes)")
    print(f"Common genes: {len(common_genes)} ({len(common_genes)/total_genes*100:.1f}% of total genes)")

def plot_top_genes_heatmap(data, n_top=50):
    """
    Create heatmap showing peak coverage for top genes.
    """
    for cell_type in data:
        for condition in data[cell_type]:
            if not data[cell_type][condition].empty:
                df = data[cell_type][condition]
                
                # Get top genes by total coverage
                top_genes = df.nlargest(n_top, 'total_coverage')
                
                # Create a normalized coverage score
                coverage_scores = top_genes[['num_peaks', 'total_coverage']].values
                normalized_scores = (coverage_scores - coverage_scores.min(axis=0)) / \
                                 (coverage_scores.max(axis=0) - coverage_scores.min(axis=0))
                
                plt.figure(figsize=(10, 12))
                sns.heatmap(normalized_scores,
                          cmap='YlOrRd',
                          xticklabels=['Peak Count', 'Coverage'],
                          yticklabels=top_genes['gene_name'],
                          cbar_kws={'label': 'Normalized Score'})
                
                plt.title(f'Top {n_top} Genes - {cell_type} {condition}')
                plt.tight_layout()
                plt.show()

def create_summary_statistics(data):
    """
    Generate and display summary statistics.
    """
    summary_stats = []
    
    for cell_type in data:
        for condition in data[cell_type]:
            if not data[cell_type][condition].empty:
                df = data[cell_type][condition]
                stats = {
                    'Cell Type': cell_type,
                    'Condition': condition,
                    'Total Genes': len(df),
                    'Mean Peaks per Gene': df['num_peaks'].mean(),
                    'Median Peaks per Gene': df['num_peaks'].median(),
                    'Mean Distance to TSS': df['min_distance_to_tss'].mean(),
                    'Median Distance to TSS': df['min_distance_to_tss'].median(),
                    'Total Peak Coverage': df['total_coverage'].sum(),
                    'Mean Peak Coverage': df['total_coverage'].mean()
                }
                summary_stats.append(stats)
    
    if summary_stats:
        stats_df = pd.DataFrame(summary_stats)
        display(stats_df)
        
        # Create summary plot
        plt.figure(figsize=(10, 6))
        sns.barplot(data=stats_df, x='Cell Type', y='Total Genes', hue='Condition')
        plt.title('Total Genes with CpG-overlapping Peaks')
        plt.tight_layout()
        plt.show()

def plot_peak_size_distribution(data):
    """
    Create boxplots showing the distribution of peak sizes.
    """
    plt.figure(figsize=(10, 6))
    peak_sizes = []
    
    for cell_type in data:
        for condition in data[cell_type]:
            if not data[cell_type][condition].empty:
                sizes = data[cell_type][condition]['avg_peak_size']
                labels = [f"{cell_type}-{condition}"] * len(sizes)
                peak_sizes.extend(list(zip(sizes, labels)))
    
    if peak_sizes:
        df = pd.DataFrame(peak_sizes, columns=['Size', 'Group'])
        sns.boxplot(data=df, x='Group', y='Size')
        plt.xticks(rotation=45)
        plt.title('Distribution of Peak Sizes')
        plt.ylabel('Peak Size (bp)')
        plt.tight_layout()
        plt.show()

def generate_all_visualizations(results_dir):
    """
    Generate and display all visualization plots.
    
    Args:
        results_dir (str): Directory containing the analysis results
    """
    # Load all data
    data = load_data(results_dir)
    
    # Generate and display all plots
    plot_peak_distribution(data)
    plot_distance_distribution(data)
    create_venn_diagrams(data)
    plot_top_genes_heatmap(data)
    create_summary_statistics(data)
    plot_peak_size_distribution(data)
    
    print("All visualizations have been displayed")

def compare_exo_endo_coverage(data):
    """
    Compare CpG islands coverage between Exo and Endo MeCP2 within each cell type.
    """
    # Compare in NSCs
    nsc_exo = data['NSC']['Exo']
    nsc_endo = data['NSC']['Endo']
    
    # Create comparison plots for NSCs
    plt.figure(figsize=(12, 6))
    
    # Plot 1: Coverage Distribution for NSCs
    plt.subplot(1, 2, 1)
    sns.kdeplot(data=nsc_exo['mean_cpg_coverage'], label='Exo', color='blue')
    sns.kdeplot(data=nsc_endo['mean_cpg_coverage'], label='Endo', color='red')
    plt.xlabel('Mean CpG Coverage (%)')
    plt.ylabel('Density')
    plt.title('NSC MeCP2 CpG Coverage Distribution')
    plt.legend()
    
    # Plot 2: Box Plot for NSCs
    plt.subplot(1, 2, 2)
    plot_data = pd.concat([
        pd.DataFrame({'Type': 'Exo', 'Coverage': nsc_exo['mean_cpg_coverage']}),
        pd.DataFrame({'Type': 'Endo', 'Coverage': nsc_endo['mean_cpg_coverage']})
    ])
    sns.boxplot(data=plot_data, x='Type', y='Coverage')
    plt.title('NSC CpG Coverage Comparison')
    
    plt.tight_layout()
    plt.show()
    
    # Statistical comparison for NSCs
    stat_result_nsc = stats.mannwhitneyu(nsc_exo['mean_cpg_coverage'],
                                        nsc_endo['mean_cpg_coverage'],
                                        alternative='two-sided')
    
    # Display statistical results for NSCs
    print("NSC Exo vs Endo Mann-Whitney U test results:")
    print(f"Statistic: {stat_result_nsc.statistic}")
    print(f"p-value: {stat_result_nsc.pvalue}")
    print("\nNSC Summary Statistics:")
    print(f"Exo mean coverage: {nsc_exo['mean_cpg_coverage'].mean():.2f}%")
    print(f"Endo mean coverage: {nsc_endo['mean_cpg_coverage'].mean():.2f}%")
    print(f"Exo median coverage: {nsc_exo['mean_cpg_coverage'].median():.2f}%")
    print(f"Endo median coverage: {nsc_endo['mean_cpg_coverage'].median():.2f}%")

    # Compare in Neurons
    neuron_exo = data['Neuron']['Exo']
    neuron_endo = data['Neuron']['Endo']
    
    # Create comparison plots for Neurons
    plt.figure(figsize=(12, 6))
    
    # Plot 1: Coverage Distribution for Neurons
    plt.subplot(1, 2, 1)
    sns.kdeplot(data=neuron_exo['mean_cpg_coverage'], label='Exo', color='blue')
    sns.kdeplot(data=neuron_endo['mean_cpg_coverage'], label='Endo', color='red')
    plt.xlabel('Mean CpG Coverage (%)')
    plt.ylabel('Density')
    plt.title('Neuron MeCP2 CpG Coverage Distribution')
    plt.legend()
    
    # Plot 2: Box Plot for Neurons
    plt.subplot(1, 2, 2)
    plot_data = pd.concat([
        pd.DataFrame({'Type': 'Exo', 'Coverage': neuron_exo['mean_cpg_coverage']}),
        pd.DataFrame({'Type': 'Endo', 'Coverage': neuron_endo['mean_cpg_coverage']})
    ])
    sns.boxplot(data=plot_data, x='Type', y='Coverage')
    plt.title('Neuron CpG Coverage Comparison')
    
    plt.tight_layout()
    plt.show()
    
    # Statistical comparison for Neurons
    stat_result_neuron = stats.mannwhitneyu(neuron_exo['mean_cpg_coverage'],
                                          neuron_endo['mean_cpg_coverage'],
                                          alternative='two-sided')
    
    # Display statistical results for Neurons
    print("\nNeuron Exo vs Endo Mann-Whitney U test results:")
    print(f"Statistic: {stat_result_neuron.statistic}")
    print(f"p-value: {stat_result_neuron.pvalue}")
    print("\nNeuron Summary Statistics:")
    print(f"Exo mean coverage: {neuron_exo['mean_cpg_coverage'].mean():.2f}%")
    print(f"Endo mean coverage: {neuron_endo['mean_cpg_coverage'].mean():.2f}%")
    print(f"Exo median coverage: {neuron_exo['mean_cpg_coverage'].median():.2f}%")
    print(f"Endo median coverage: {neuron_endo['mean_cpg_coverage'].median():.2f}%")

########################################################################################################################33

def compare_endo_cpg_coverage_common(data):
    """
    Compare CpG islands coverage of Endo MeCP2 between NSCs and Neurons,
    but only for genes common to both cell types.
    """
    # Extract Endo data for both cell types
    nsc_endo = data['NSC']['Endo']
    neuron_endo = data['Neuron']['Endo']
    
    # Find common genes
    common_genes = set(nsc_endo['gene_name']).intersection(set(neuron_endo['gene_name']))
    print(f"Number of common genes: {len(common_genes)}")
    
    # Filter data for common genes
    nsc_endo_common = nsc_endo[nsc_endo['gene_name'].isin(common_genes)]
    neuron_endo_common = neuron_endo[neuron_endo['gene_name'].isin(common_genes)]
    
    # Create comparison plots
    plt.figure(figsize=(12, 6))
    
    # Plot 1: Coverage Distribution
    plt.subplot(1, 2, 1)
    sns.kdeplot(data=nsc_endo_common['mean_cpg_coverage'], label='NSCs', color='#bcdae6', linewidth=3.0)
    sns.kdeplot(data=neuron_endo_common['mean_cpg_coverage'], label='Neurons', color='#efb3b1', linewidth=3.0)
    plt.xlabel('Mean CpG Coverage (%)')
    plt.ylabel('Density')
    plt.title('Endo MeCP2 CpG Coverage Distribution\n(Common Genes)')
    plt.legend()
    
    # Plot 2: Box Plot
    plt.subplot(1, 2, 2)
    plot_data = pd.concat([
        pd.DataFrame({'Cell Type': 'NSCs', 'Coverage': nsc_endo_common['mean_cpg_coverage']}),
        pd.DataFrame({'Cell Type': 'Neurons', 'Coverage': neuron_endo_common['mean_cpg_coverage']})
    ])
    sns.boxplot(data=plot_data, x='Cell Type', y='Coverage')
    plt.title('CpG Coverage Comparison\n(Common Genes)')
    
    plt.tight_layout()
    plt.show()
    
    # Statistical comparison
    stat_result = stats.mannwhitneyu(nsc_endo_common['mean_cpg_coverage'],
                                   neuron_endo_common['mean_cpg_coverage'],
                                   alternative='two-sided')
    
    # Display statistical results
    print("\nMann-Whitney U test results (Common Genes):")
    print(f"Statistic: {stat_result.statistic}")
    print(f"p-value: {stat_result.pvalue}")
    print("\nSummary Statistics (Common Genes):")
    print(f"NSCs mean coverage: {nsc_endo_common['mean_cpg_coverage'].mean():.2f}%")
    print(f"Neurons mean coverage: {neuron_endo_common['mean_cpg_coverage'].mean():.2f}%")
    print(f"NSCs median coverage: {nsc_endo_common['mean_cpg_coverage'].median():.2f}%")
    print(f"Neurons median coverage: {neuron_endo_common['mean_cpg_coverage'].median():.2f}%")

    # Do the same for Exo data
    nsc_exo = data['NSC']['Exo']
    neuron_exo = data['Neuron']['Exo']
    
    # Find common genes for Exo
    common_genes_exo = set(nsc_exo['gene_name']).intersection(set(neuron_exo['gene_name']))
    print(f"\nNumber of common genes (Exo): {len(common_genes_exo)}")
    
    # Filter data for common genes
    nsc_exo_common = nsc_exo[nsc_exo['gene_name'].isin(common_genes_exo)]
    neuron_exo_common = neuron_exo[neuron_exo['gene_name'].isin(common_genes_exo)]
    
    plt.figure(figsize=(12, 6))
    
    # Plot 1: Coverage Distribution for Exo
    plt.subplot(1, 2, 1)
    sns.kdeplot(data=nsc_exo_common['mean_cpg_coverage'], label='NSCs', color='#1919f5', linewidth=3.0)
    sns.kdeplot(data=neuron_exo_common['mean_cpg_coverage'], label='Neurons', color='#ec4b3d', linewidth=3.0)
    plt.xlabel('Mean CpG Coverage (%)')
    plt.ylabel('Density')
    plt.title('Exo MeCP2 CpG Coverage Distribution\n(Common Genes)')
    plt.legend()
    
    # Plot 2: Box Plot for Exo
    plt.subplot(1, 2, 2)
    plot_data = pd.concat([
        pd.DataFrame({'Cell Type': 'NSCs', 'Coverage': nsc_exo_common['mean_cpg_coverage']}),
        pd.DataFrame({'Cell Type': 'Neurons', 'Coverage': neuron_exo_common['mean_cpg_coverage']})
    ])
    sns.boxplot(data=plot_data, x='Cell Type', y='Coverage')
    plt.title('CpG Coverage Comparison\n(Common Genes)')
    
    plt.tight_layout()
    plt.show()
    
    # Statistical comparison for Exo
    stat_result_exo = stats.mannwhitneyu(nsc_exo_common['mean_cpg_coverage'],
                                       neuron_exo_common['mean_cpg_coverage'],
                                       alternative='two-sided')
    
    # Display statistical results for Exo
    print("\nExo MeCP2 Mann-Whitney U test results (Common Genes):")
    print(f"Statistic: {stat_result_exo.statistic}")
    print(f"p-value: {stat_result_exo.pvalue}")
    print("\nSummary Statistics (Common Genes):")
    print(f"NSCs mean coverage: {nsc_exo_common['mean_cpg_coverage'].mean():.2f}%")
    print(f"Neurons mean coverage: {neuron_exo_common['mean_cpg_coverage'].mean():.2f}%")
    print(f"NSCs median coverage: {nsc_exo_common['mean_cpg_coverage'].median():.2f}%")
    print(f"Neurons median coverage: {neuron_exo_common['mean_cpg_coverage'].median():.2f}%")


def compare_exo_endo_coverage_common(data):
    """
    Compare CpG islands coverage between Exo and Endo MeCP2 within each cell type,
    considering only genes common to both Exo and Endo conditions.
    """
    # Compare in NSCs
    nsc_exo = data['NSC']['Exo']
    nsc_endo = data['NSC']['Endo']
    
    # Find common genes in NSCs
    nsc_common_genes = set(nsc_exo['gene_name']).intersection(set(nsc_endo['gene_name']))
    print(f"Number of common genes in NSCs: {len(nsc_common_genes)}")
    
    # Filter for common genes
    nsc_exo_common = nsc_exo[nsc_exo['gene_name'].isin(nsc_common_genes)]
    nsc_endo_common = nsc_endo[nsc_endo['gene_name'].isin(nsc_common_genes)]
    
    # Create comparison plots for NSCs
    plt.figure(figsize=(12, 6))
    
    # Plot 1: Coverage Distribution for NSCs
    plt.subplot(1, 2, 1)
    sns.kdeplot(data=nsc_exo_common['mean_cpg_coverage'], label='Exo', color='#1919f5', linewidth=3.0)
    sns.kdeplot(data=nsc_endo_common['mean_cpg_coverage'], label='Endo', color='#bcdae6', linewidth=3.0)
    plt.xlabel('Mean CpG Coverage (%)')
    plt.ylabel('Density')
    plt.title('NSC MeCP2 CpG Coverage Distribution\n(Common Genes)')
    plt.legend()
    
    # Plot 2: Box Plot for NSCs
    plt.subplot(1, 2, 2)
    plot_data = pd.concat([
        pd.DataFrame({'Type': 'Exo', 'Coverage': nsc_exo_common['mean_cpg_coverage']}),
        pd.DataFrame({'Type': 'Endo', 'Coverage': nsc_endo_common['mean_cpg_coverage']})
    ])
    sns.boxplot(data=plot_data, x='Type', y='Coverage')
    plt.title('NSC CpG Coverage Comparison\n(Common Genes)')
    
    plt.tight_layout()
    plt.show()
    
    # Statistical comparison for NSCs
    stat_result_nsc = stats.mannwhitneyu(nsc_exo_common['mean_cpg_coverage'],
                                        nsc_endo_common['mean_cpg_coverage'],
                                        alternative='two-sided')
    
    # Display statistical results for NSCs
    print("\nNSC Exo vs Endo Mann-Whitney U test results (Common Genes):")
    print(f"Statistic: {stat_result_nsc.statistic}")
    print(f"p-value: {stat_result_nsc.pvalue}")
    print("\nNSC Summary Statistics (Common Genes):")
    print(f"Exo mean coverage: {nsc_exo_common['mean_cpg_coverage'].mean():.2f}%")
    print(f"Endo mean coverage: {nsc_endo_common['mean_cpg_coverage'].mean():.2f}%")
    print(f"Exo median coverage: {nsc_exo_common['mean_cpg_coverage'].median():.2f}%")
    print(f"Endo median coverage: {nsc_endo_common['mean_cpg_coverage'].median():.2f}%")

    # Compare in Neurons
    neuron_exo = data['Neuron']['Exo']
    neuron_endo = data['Neuron']['Endo']
    
    # Find common genes in Neurons
    neuron_common_genes = set(neuron_exo['gene_name']).intersection(set(neuron_endo['gene_name']))
    print(f"\nNumber of common genes in Neurons: {len(neuron_common_genes)}")
    
    # Filter for common genes
    neuron_exo_common = neuron_exo[neuron_exo['gene_name'].isin(neuron_common_genes)]
    neuron_endo_common = neuron_endo[neuron_endo['gene_name'].isin(neuron_common_genes)]
    
    # Create comparison plots for Neurons
    plt.figure(figsize=(12, 6))
    
    # Plot 1: Coverage Distribution for Neurons
    plt.subplot(1, 2, 1)
    sns.kdeplot(data=neuron_exo_common['mean_cpg_coverage'], label='Exo', color='#ec4b3d', linewidth=3.0)
    sns.kdeplot(data=neuron_endo_common['mean_cpg_coverage'], label='Endo', color='#efb3b1', linewidth=3.0)
    plt.xlabel('Mean CpG Coverage (%)')
    plt.ylabel('Density')
    plt.title('Neuron MeCP2 CpG Coverage Distribution\n(Common Genes)')
    plt.legend()
    
    # Plot 2: Box Plot for Neurons
    plt.subplot(1, 2, 2)
    plot_data = pd.concat([
        pd.DataFrame({'Type': 'Exo', 'Coverage': neuron_exo_common['mean_cpg_coverage']}),
        pd.DataFrame({'Type': 'Endo', 'Coverage': neuron_endo_common['mean_cpg_coverage']})
    ])
    sns.boxplot(data=plot_data, x='Type', y='Coverage')
    plt.title('Neuron CpG Coverage Comparison\n(Common Genes)')
    
    plt.tight_layout()
    plt.show()
    
    # Statistical comparison for Neurons
    stat_result_neuron = stats.mannwhitneyu(neuron_exo_common['mean_cpg_coverage'],
                                          neuron_endo_common['mean_cpg_coverage'],
                                          alternative='two-sided')
    
    # Display statistical results for Neurons
    print("\nNeuron Exo vs Endo Mann-Whitney U test results (Common Genes):")
    print(f"Statistic: {stat_result_neuron.statistic}")
    print(f"p-value: {stat_result_neuron.pvalue}")
    print("\nNeuron Summary Statistics (Common Genes):")
    print(f"Exo mean coverage: {neuron_exo_common['mean_cpg_coverage'].mean():.2f}%")
    print(f"Endo mean coverage: {neuron_endo_common['mean_cpg_coverage'].mean():.2f}%")
    print(f"Exo median coverage: {neuron_exo_common['mean_cpg_coverage'].median():.2f}%")
    print(f"Endo median coverage: {neuron_endo_common['mean_cpg_coverage'].median():.2f}%")



########################################### Legacy code ###########################################

# %%script false --no-raise-error

# def load_and_process_data(results_dir):
#     """
#     Load and process data with additional CpG coverage calculations.
#     """
#     cell_types = ['NSC', 'Neuron']
#     conditions = ['Exo', 'Endo']
#     data = {}
    
#     for cell_type in cell_types:
#         data[cell_type] = {}
#         for condition in conditions:
#             file_path = os.path.join(results_dir, f"{cell_type}_{condition}_cpg_genes.tsv")
#             if os.path.exists(file_path):
#                 df = pd.read_csv(file_path, sep='\t')
#                 # Calculate CpG coverage percentage
#                 if 'total_coverage' in df.columns:
#                     df['cpg_coverage_percent'] = (df['total_coverage'] / df['total_coverage'].max()) * 100
#                 data[cell_type][condition] = df
#             else:
#                 print(f"Warning: Missing file {file_path}")
#                 data[cell_type][condition] = pd.DataFrame()
    
#     return data

# def compare_endo_cpg_coverage(data):
#     """
#     Compare CpG islands coverage of Endo MeCP2 between NPCs and Neurons.
#     """
#     # Extract Endo data for both cell types
#     npc_endo = data['NSC']['Endo']
#     neuron_endo = data['Neuron']['Endo']
    
#     # Create comparison plots
#     plt.figure(figsize=(12, 6))
    
#     # Plot 1: Coverage Distribution
#     plt.subplot(1, 2, 1)
#     sns.kdeplot(data=npc_endo['cpg_coverage_percent'], label='NPCs', color='blue')
#     sns.kdeplot(data=neuron_endo['cpg_coverage_percent'], label='Neurons', color='red')
#     plt.xlabel('CpG Coverage (%)')
#     plt.ylabel('Density')
#     plt.title('Endo MeCP2 CpG Coverage Distribution')
#     plt.legend()
    
#     # Plot 2: Box Plot
#     plt.subplot(1, 2, 2)
#     plot_data = pd.concat([
#         pd.DataFrame({'Cell Type': 'NPCs', 'Coverage': npc_endo['cpg_coverage_percent']}),
#         pd.DataFrame({'Cell Type': 'Neurons', 'Coverage': neuron_endo['cpg_coverage_percent']})
#     ])
#     sns.boxplot(data=plot_data, x='Cell Type', y='Coverage')
#     plt.title('CpG Coverage Comparison')
    
#     plt.tight_layout()
#     plt.show()
    
#     # Statistical comparison using the correct function name
#     stat_result = stats.mannwhitneyu(npc_endo['cpg_coverage_percent'],
#                                    neuron_endo['cpg_coverage_percent'],
#                                    alternative='two-sided')
    
#     # Display statistical results
#     print("Mann-Whitney U test results:")
#     print(f"Statistic: {stat_result.statistic}")
#     print(f"p-value: {stat_result.pvalue}")
#     print("\nSummary Statistics:")
#     print(f"NPCs mean coverage: {npc_endo['cpg_coverage_percent'].mean():.2f}%")
#     print(f"Neurons mean coverage: {neuron_endo['cpg_coverage_percent'].mean():.2f}%")
#     print(f"NPCs median coverage: {npc_endo['cpg_coverage_percent'].median():.2f}%")
#     print(f"Neurons median coverage: {neuron_endo['cpg_coverage_percent'].median():.2f}%")

#     # Create identical plots for Exo comparison
#     npc_exo = data['NSC']['Exo']
#     neuron_exo = data['Neuron']['Exo']
    
#     plt.figure(figsize=(12, 6))
    
#     # Plot 1: Coverage Distribution for Exo
#     plt.subplot(1, 2, 1)
#     sns.kdeplot(data=npc_exo['cpg_coverage_percent'], label='NPCs', color='blue')
#     sns.kdeplot(data=neuron_exo['cpg_coverage_percent'], label='Neurons', color='red')
#     plt.xlabel('CpG Coverage (%)')
#     plt.ylabel('Density')
#     plt.title('Exo MeCP2 CpG Coverage Distribution')
#     plt.legend()
    
#     # Plot 2: Box Plot for Exo
#     plt.subplot(1, 2, 2)
#     plot_data = pd.concat([
#         pd.DataFrame({'Cell Type': 'NPCs', 'Coverage': npc_exo['cpg_coverage_percent']}),
#         pd.DataFrame({'Cell Type': 'Neurons', 'Coverage': neuron_exo['cpg_coverage_percent']})
#     ])
#     sns.boxplot(data=plot_data, x='Cell Type', y='Coverage')
#     plt.title('CpG Coverage Comparison')
    
#     plt.tight_layout()
#     plt.show()
    
#     # Statistical comparison for Exo
#     stat_result_exo = stats.mannwhitneyu(npc_exo['cpg_coverage_percent'],
#                                        neuron_exo['cpg_coverage_percent'],
#                                        alternative='two-sided')
    
#     # Display statistical results for Exo
#     print("\nExo MeCP2 Mann-Whitney U test results:")
#     print(f"Statistic: {stat_result_exo.statistic}")
#     print(f"p-value: {stat_result_exo.pvalue}")
#     print("\nSummary Statistics:")
#     print(f"NPCs mean coverage: {npc_exo['cpg_coverage_percent'].mean():.2f}%")
#     print(f"Neurons mean coverage: {neuron_exo['cpg_coverage_percent'].mean():.2f}%")
#     print(f"NPCs median coverage: {npc_exo['cpg_coverage_percent'].median():.2f}%")
#     print(f"Neurons median coverage: {neuron_exo['cpg_coverage_percent'].median():.2f}%")

# def analyze_common_cpg_targets(data):
#     """
#     Analyze CpG targets common between all four experimental groups.
#     """
#     # Extract gene sets
#     npc_endo = set(data['NSC']['Endo']['gene_name'])
#     neuron_endo = set(data['Neuron']['Endo']['gene_name'])
#     npc_exo = set(data['NSC']['Exo']['gene_name'])
#     neuron_exo = set(data['Neuron']['Exo']['gene_name'])
    
#     # Calculate common targets across all groups
#     all_common = npc_endo.intersection(neuron_endo, npc_exo, neuron_exo)
#     all_unique = npc_endo.union(neuron_endo, npc_exo, neuron_exo)
    
#     # Calculate percentage of common targets
#     common_percent = (len(all_common) / len(all_unique)) * 100
    
#     # Create the 4-way Venn diagram
#     plt.figure(figsize=(10, 10))
#     venn({
#         'NSC Endo': npc_endo,
#         'Neuron Endo': neuron_endo,
#         'NSC Exo': npc_exo,
#         'Neuron Exo': neuron_exo
#     })
#     plt.title('CpG Targets Across All Groups')
#     plt.tight_layout()
#     plt.show()
    
#     # Display overlap statistics
#     print("CpG Target Statistics:")
#     print("-" * 50)
#     print(f"NSC Endo targets: {len(npc_endo)}")
#     print(f"Neuron Endo targets: {len(neuron_endo)}")
#     print(f"NSC Exo targets: {len(npc_exo)}")
#     print(f"Neuron Exo targets: {len(neuron_exo)}")
#     print("-" * 50)
#     print(f"Total unique targets: {len(all_unique)}")
#     print(f"Common targets (all groups): {len(all_common)} ({common_percent:.2f}%)")
    
#     # Calculate pairwise overlaps
#     pairs = [
#         ('NSC Endo', 'Neuron Endo', npc_endo, neuron_endo),
#         ('NSC Endo', 'NSC Exo', npc_endo, npc_exo),
#         ('NSC Endo', 'Neuron Exo', npc_endo, neuron_exo),
#         ('Neuron Endo', 'NSC Exo', neuron_endo, npc_exo),
#         ('Neuron Endo', 'Neuron Exo', neuron_endo, neuron_exo),
#         ('NSC Exo', 'Neuron Exo', npc_exo, neuron_exo)
#     ]
    
#     print("\nPairwise Overlaps:")
#     print("-" * 50)
#     for label1, label2, set1, set2 in pairs:
#         overlap = len(set1.intersection(set2))
#         total = len(set1.union(set2))
#         percent = (overlap / total) * 100
#         print(f"{label1} vs {label2}: {overlap} ({percent:.2f}%)")
    
#     # Calculate group-specific targets
#     print("\nGroup-Specific Targets:")
#     print("-" * 50)
#     nsc_endo_specific = len(npc_endo - (neuron_endo | npc_exo | neuron_exo))
#     neuron_endo_specific = len(neuron_endo - (npc_endo | npc_exo | neuron_exo))
#     nsc_exo_specific = len(npc_exo - (npc_endo | neuron_endo | neuron_exo))
#     neuron_exo_specific = len(neuron_exo - (npc_endo | neuron_endo | npc_exo))
    
#     print(f"NSC Endo-specific: {nsc_endo_specific}")
#     print(f"Neuron Endo-specific: {neuron_endo_specific}")
#     print(f"NSC Exo-specific: {nsc_exo_specific}")
#     print(f"Neuron Exo-specific: {neuron_exo_specific}")

# def analyze_exo_vs_endo_enrichment(data):
#     """
#     Analyze enrichment of Exo vs Endo MeCP2 CpG coverage.
#     """
#     enriched_genes = {'NSC': [], 'Neuron': []}
    
#     for cell_type in ['NSC', 'Neuron']:
#         exo_df = data[cell_type]['Exo']
#         endo_df = data[cell_type]['Endo']
        
#         # Merge Exo and Endo data
#         merged = pd.merge(exo_df, endo_df, 
#                          on='gene_name', 
#                          suffixes=('_exo', '_endo'))
        
#         # Calculate enrichment ratio
#         merged['enrichment_ratio'] = merged['total_coverage_exo'] / merged['total_coverage_endo']
        
#         # Filter for significantly enriched genes (log2 fold change >= 0.5)
#         merged['log2_enrichment'] = np.log2(merged['enrichment_ratio'])
#         enriched = merged[merged['log2_enrichment'] >= 0.5].sort_values('enrichment_ratio', ascending=False)
#         enriched_genes[cell_type] = enriched
        
#         # Save enriched genes to TSV files
#         output_file = f"results/enriched_genes_{cell_type}.tsv"
#         enriched.to_csv(output_file, sep='\t', index=False)
        
#         # Display enriched genes
#         # print(f"\nTop enriched genes in {cell_type}:")
#         # display(enriched.head())
    
#     # Find common enriched genes
#     common_enriched = set(enriched_genes['NSC']['gene_name']).intersection(
#         set(enriched_genes['Neuron']['gene_name']))
    
#     print("\nGenes enriched in Exo vs Endo in both NPCs and Neurons:")
#     for gene in sorted(common_enriched):
#         print(gene)

# def generate_comprehensive_analysis(results_dir):
#     """
#     Generate all analyses and visualizations.
#     """
#     # Load data
#     data = load_and_process_data(results_dir)
    
#     # Generate all analyses
#     compare_endo_cpg_coverage(data)
#     analyze_common_cpg_targets(data)
#     analyze_exo_vs_endo_enrichment(data)


# %%script false --no-raise-error

# def load_data(results_dir):
#     """
#     Load all result files into a dictionary of DataFrames.
    
#     Args:
#         results_dir (str): Path to the results directory
    
#     Returns:
#         dict: Nested dictionary containing DataFrames for each cell type and condition
#     """
#     cell_types = ['NSC', 'Neuron']
#     conditions = ['Exo', 'Endo']
#     data = {}
    
#     for cell_type in cell_types:
#         data[cell_type] = {}
#         for condition in conditions:
#             file_path = os.path.join(results_dir, f"{cell_type}_{condition}_cpg_genes.tsv")
#             if os.path.exists(file_path):
#                 data[cell_type][condition] = pd.read_csv(file_path, sep='\t')
#             else:
#                 print(f"Warning: Missing file {file_path}")
#                 data[cell_type][condition] = pd.DataFrame()
                
#     return data

# def plot_peak_distribution(data):
#     """
#     Create violin plots showing the distribution of peaks per gene.
#     """
#     plt.figure(figsize=(10, 6))
#     data_points = []
    
#     for cell_type in data:
#         for condition in data[cell_type]:
#             if not data[cell_type][condition].empty:
#                 peak_counts = data[cell_type][condition]['num_peaks']
#                 label = f"{cell_type}-{condition}"
#                 data_points.extend([(x, label) for x in peak_counts])
    
#     if data_points:
#         df = pd.DataFrame(data_points, columns=['Peaks', 'Group'])
#         sns.violinplot(data=df, x='Group', y='Peaks')
#         plt.xticks(rotation=45)
#         plt.title('Distribution of CpG-overlapping Peaks per Gene')
#         plt.tight_layout()
#         plt.show()

# def plot_distance_distribution(data):
#     """
#     Create KDE plots showing the distribution of distances to TSS.
#     """
#     plt.figure(figsize=(12, 6))
#     colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
#     color_idx = 0
    
#     for cell_type in data:
#         for condition in data[cell_type]:
#             if not data[cell_type][condition].empty:
#                 distances = data[cell_type][condition]['min_distance_to_tss']
#                 sns.kdeplot(data=distances, label=f"{cell_type}-{condition}",
#                           color=colors[color_idx])
#                 color_idx += 1
    
#     plt.xlabel('Distance to TSS (bp)')
#     plt.ylabel('Density')
#     plt.title('Distribution of Minimum Distances to TSS')
#     plt.legend()
#     plt.tight_layout()
#     plt.show()

# def create_venn_diagrams(data):
#     """
#     Create Venn diagrams showing overlap between conditions for each cell type.
#     """
#     colors = {'NSC': ['#1f77b4', '#2ca02c'], 'Neuron': ['#ff7f0e', '#d62728']}
    
#     for cell_type in data:
#         plt.figure(figsize=(8, 8))
        
#         exo_genes = set(data[cell_type]['Exo']['gene_name']) if not data[cell_type]['Exo'].empty else set()
#         endo_genes = set(data[cell_type]['Endo']['gene_name']) if not data[cell_type]['Endo'].empty else set()
        
#         # Create Venn diagram
#         venn2([exo_genes, endo_genes],
#               set_labels=('Exogenous', 'Endogenous'),
#               set_colors=colors[cell_type])
        
#         # Add title separately
#         plt.title(f'{cell_type} Gene Overlap')
        
#         plt.tight_layout()
#         plt.show()

# def plot_top_genes_heatmap(data, n_top=50):
#     """
#     Create heatmap showing peak coverage for top genes.
#     """
#     for cell_type in data:
#         for condition in data[cell_type]:
#             if not data[cell_type][condition].empty:
#                 df = data[cell_type][condition]
                
#                 # Get top genes by total coverage
#                 top_genes = df.nlargest(n_top, 'total_coverage')
                
#                 # Create a normalized coverage score
#                 coverage_scores = top_genes[['num_peaks', 'total_coverage']].values
#                 normalized_scores = (coverage_scores - coverage_scores.min(axis=0)) / \
#                                  (coverage_scores.max(axis=0) - coverage_scores.min(axis=0))
                
#                 plt.figure(figsize=(10, 12))
#                 sns.heatmap(normalized_scores,
#                           cmap='YlOrRd',
#                           xticklabels=['Peak Count', 'Coverage'],
#                           yticklabels=top_genes['gene_name'],
#                           cbar_kws={'label': 'Normalized Score'})
                
#                 plt.title(f'Top {n_top} Genes - {cell_type} {condition}')
#                 plt.tight_layout()
#                 plt.show()

# def create_summary_statistics(data):
#     """
#     Generate and display summary statistics.
#     """
#     summary_stats = []
    
#     for cell_type in data:
#         for condition in data[cell_type]:
#             if not data[cell_type][condition].empty:
#                 df = data[cell_type][condition]
#                 stats = {
#                     'Cell Type': cell_type,
#                     'Condition': condition,
#                     'Total Genes': len(df),
#                     'Mean Peaks per Gene': df['num_peaks'].mean(),
#                     'Median Peaks per Gene': df['num_peaks'].median(),
#                     'Mean Distance to TSS': df['min_distance_to_tss'].mean(),
#                     'Median Distance to TSS': df['min_distance_to_tss'].median(),
#                     'Total Peak Coverage': df['total_coverage'].sum(),
#                     'Mean Peak Coverage': df['total_coverage'].mean()
#                 }
#                 summary_stats.append(stats)
    
#     if summary_stats:
#         stats_df = pd.DataFrame(summary_stats)
#         display(stats_df)
        
#         # Create summary plot
#         plt.figure(figsize=(10, 6))
#         sns.barplot(data=stats_df, x='Cell Type', y='Total Genes', hue='Condition')
#         plt.title('Total Genes with CpG-overlapping Peaks')
#         plt.tight_layout()
#         plt.show()

# def plot_peak_size_distribution(data):
#     """
#     Create boxplots showing the distribution of peak sizes.
#     """
#     plt.figure(figsize=(10, 6))
#     peak_sizes = []
    
#     for cell_type in data:
#         for condition in data[cell_type]:
#             if not data[cell_type][condition].empty:
#                 sizes = data[cell_type][condition]['avg_peak_size']
#                 labels = [f"{cell_type}-{condition}"] * len(sizes)
#                 peak_sizes.extend(list(zip(sizes, labels)))
    
#     if peak_sizes:
#         df = pd.DataFrame(peak_sizes, columns=['Size', 'Group'])
#         sns.boxplot(data=df, x='Group', y='Size')
#         plt.xticks(rotation=45)
#         plt.title('Distribution of Peak Sizes')
#         plt.ylabel('Peak Size (bp)')
#         plt.tight_layout()
#         plt.show()

# def generate_all_visualizations(results_dir):
#     """
#     Generate and display all visualization plots.
    
#     Args:
#         results_dir (str): Directory containing the analysis results
#     """
#     # Load all data
#     data = load_data(results_dir)
    
#     # Generate and display all plots
#     plot_peak_distribution(data)
#     plot_distance_distribution(data)
#     create_venn_diagrams(data)
#     plot_top_genes_heatmap(data)
#     create_summary_statistics(data)
#     plot_peak_size_distribution(data)
    
#     print("All visualizations have been displayed")



# %%script false --no-raise-error
# 
# def get_peaks_with_cpg(peak_file, cpg_file, extend=300, coverage_threshold=20):
#     """
#     Identifies peaks that overlap with CpG islands above a coverage threshold.
    
#     Args:
#         peak_file: BED file containing peak regions
#         cpg_file: BED file containing CpG islands
#         extend: Number of base pairs to extend CpG islands (default: 300bp)
#         coverage_threshold: Minimum percentage overlap required (default: 20%)
    
#     Returns:
#         set: Peak coordinates that meet the CpG overlap criteria
#     """
#     print(f"\nProcessing {peak_file}")
    
#     # Validate input files
#     if not all(os.path.exists(f) for f in [peak_file, cpg_file]):
#         print("Error: Input file(s) missing!")
#         return set()
    
#     # Filter for standard chromosomes (chr1-22, X, Y)
#     standard_chroms = [f"chr{i}" for i in range(1, 23)] + ['chrX', 'chrY']
    
#     # Create temporary filtered files
#     filtered_peaks = "temp_filtered_peaks.bed"
#     filtered_cpg = "temp_filtered_cpg.bed"
    
#     # Filter peaks and count
#     with open(peak_file) as fin, open(filtered_peaks, 'w') as fout:
#         peaks = [line for line in fin if line.split('\t')[0] in standard_chroms]
#         fout.writelines(peaks)
#         print(f"Peaks after chromosome filtering: {len(peaks)}")
    
#     # Filter and extend CpG islands
#     with open(cpg_file) as fin, open(filtered_cpg, 'w') as fout:
#         fout.writelines(line for line in fin if line.split('\t')[0] in standard_chroms)
    
#     # # Print first 5 lines from each file for debugging
#     # with open(filtered_peaks) as f:
#     #     print("\nFirst 5 lines of filtered peaks:")
#     #     for i, line in enumerate(f):
#     #         if i < 5:
#     #             print(line.strip())  
#     # with open(filtered_cpg) as f:
#     #     print("\nFirst 5 lines of filtered CpG:")
#     #     for i, line in enumerate(f):
#     #         if i < 5:
#     #             print(line.strip())
    
#     # Extend CpG regions and find overlaps using bedtools
#     subprocess.run(f"bedtools slop -i {filtered_cpg} -g DATA/genome.size -b {extend} > temp_extended_cpg.bed", shell=True)
#     result = subprocess.run(f"bedtools intersect -a {filtered_peaks} -b temp_extended_cpg.bed -wao", 
#                           shell=True, capture_output=True, text=True)
    
#     # # Debug: Print first few lines of bedtools output
#     # print("\nFirst 5 lines of bedtools intersect output:")
#     # for i, line in enumerate(result.stdout.strip().split('\n')[:5]):
#     #     print(line)
    
#     # Process overlaps and calculate coverage
#     peaks_with_cpg = set()
#     coverage_stats = []
#     current_peak = {'id': None, 'overlaps': []}
#     qualified_peaks = 0
#     total_peaks_with_overlap = 0
    
#     for line in result.stdout.strip().split('\n'):
#         fields = line.split('\t')
        
#         # Fields from peak (first BED file)
#         peak_chr = fields[0]
#         peak_start = fields[1]
#         peak_end = fields[2]
#         peak_id = f"{peak_chr}:{peak_start}-{peak_end}"
        
#         # Check if there's a valid CpG overlap
#         if len(fields) > 6 and fields[4] != "." and fields[4] != "-1":  # Changed condition
#             try:
#                 cpg_start = int(fields[5])  # Changed from fields[6]
#                 cpg_end = int(fields[6])    # Changed from fields[7]
#                 overlap_amount = int(fields[-1])  # Last field is overlap amount
                
#                 # If we encounter a new peak
#                 if peak_id != current_peak['id']:
#                     # Process previous peak if it exists
#                     if current_peak['overlaps']:
#                         max_coverage = max(current_peak['overlaps'])
#                         total_peaks_with_overlap += 1
#                         if max_coverage >= coverage_threshold:
#                             coverage_stats.append(max_coverage)
#                             peaks_with_cpg.add(current_peak['id'])
#                             qualified_peaks += 1
                    
#                     current_peak = {'id': peak_id, 'overlaps': []}
                
#                 # Calculate overlap coverage using CpG length
#                 if overlap_amount > 0:  # If there is overlap
#                     cpg_length = cpg_end - cpg_start  # CpG island length
#                     if cpg_length > 0:  # Avoid division by zero
#                         coverage = (overlap_amount / cpg_length) * 100
#                         current_peak['overlaps'].append(coverage)
                        
#                         # # Debug: Print coverage calculation for first few peaks
#                         # if qualified_peaks < 3:
#                         #     print(f"\nCoverage calculation for {peak_id}:")
#                         #     print(f"CpG region: {cpg_start}-{cpg_end}")
#                         #     print(f"Overlap amount: {overlap_amount}")
#                         #     print(f"CpG length: {cpg_length}")
#                         #     print(f"Coverage: {coverage:.2f}%")
            
#             except (ValueError, IndexError) as e:
#                 continue
        
#         else:
#             # If we encounter a new peak with no overlap
#             if peak_id != current_peak['id']:
#                 # Process previous peak if it exists
#                 if current_peak['overlaps']:
#                     max_coverage = max(current_peak['overlaps'])
#                     total_peaks_with_overlap += 1
#                     if max_coverage >= coverage_threshold:
#                         coverage_stats.append(max_coverage)
#                         peaks_with_cpg.add(current_peak['id'])
#                         qualified_peaks += 1
                
#                 current_peak = {'id': peak_id, 'overlaps': []}
    
#     # Process the last peak
#     if current_peak['overlaps']:
#         max_coverage = max(current_peak['overlaps'])
#         total_peaks_with_overlap += 1
#         if max_coverage >= coverage_threshold:
#             coverage_stats.append(max_coverage)
#             peaks_with_cpg.add(current_peak['id'])
#             qualified_peaks += 1
    
#     # Print coverage statistics
#     if coverage_stats:
#         print(f"\nPeaks with any CpG overlap: {total_peaks_with_overlap}")
#         print(f"Peaks meeting {coverage_threshold}% coverage threshold: {qualified_peaks}")
#         print(f"Coverage range: {min(coverage_stats):.2f}% - {max(coverage_stats):.2f}%")
#         print(f"Mean coverage of qualified peaks: {sum(coverage_stats)/len(coverage_stats):.2f}%")
#     else:
#         print("\nNo peaks met the coverage criteria")
    
#     # Cleanup temporary files
#     for f in [filtered_peaks, filtered_cpg, "temp_extended_cpg.bed"]:
#         if os.path.exists(f):
#             os.remove(f)
    
#     return peaks_with_cpg

# def get_genes_with_cpg_enrichment(peak_file, cpg_file, gtf_file, output_dir, cell_type, condition, 
#                                  extend_cpg=300, extend_tss=2000, coverage_threshold=20):
#     """
#     Identifies genes with CpG-overlapping peaks near their TSS regions.
#     """
#     print(f"\nAnalyzing {cell_type} {condition} peaks")
    
#     # Validate input files
#     if not os.path.exists(gtf_file):
#         print(f"Error: GTF file not found: {gtf_file}")
#         return None
        
#     os.makedirs(output_dir, exist_ok=True)
    
#     # 1. Get peaks overlapping with CpG islands
#     cpg_peaks = get_peaks_with_cpg(peak_file, cpg_file, extend=extend_cpg, 
#                                   coverage_threshold=coverage_threshold)
    
#     if not cpg_peaks:
#         print(f"No CpG-overlapping peaks found for {cell_type} {condition}")
#         return {'genes': pd.DataFrame(), 'total_peaks': 0}
    
#     # 2. Extract and extend TSS regions
#     try:
#         tss_bed = os.path.join(output_dir, "temp_tss.bed")
#         extract_tss_regions(gtf_file, tss_bed)
        
#         extended_tss = os.path.join(output_dir, "temp_extended_tss.bed")
#         subprocess.run(f"bedtools slop -i {tss_bed} -g DATA/genome.size -b {extend_tss} > {extended_tss}",
#                       shell=True, check=True)
        
#         # 3. Write CpG-overlapping peaks to temporary file
#         temp_peaks = os.path.join(output_dir, "temp_cpg_peaks.bed")
#         with open(temp_peaks, 'w') as f:
#             for peak in cpg_peaks:
#                 chrom, pos = peak.split(':')
#                 start, end = pos.split('-')
#                 f.write(f"{chrom}\t{start}\t{end}\n")
        
#         # 4. Find overlaps between peaks and TSS regions
#         result = subprocess.run(
#             f"bedtools intersect -a {extended_tss} -b {temp_peaks} -wo",
#             shell=True, capture_output=True, text=True, check=True
#         )
        
#         # 5. Process results
#         gene_data = defaultdict(lambda: {
#             'peaks': set(),
#             'total_peak_coverage': 0,
#             'distance_to_tss': []
#         })
        
#         for line in result.stdout.strip().split('\n'):
#             if not line:
#                 continue
            
#             fields = line.split('\t')
#             gene_name = fields[3]
#             tss_pos = int(fields[1]) + extend_tss
#             peak_start = int(fields[5])
#             peak_end = int(fields[6])
#             peak_id = f"{fields[4]}:{fields[5]}-{fields[6]}"
            
#             peak_center = (peak_start + peak_end) // 2
#             distance = abs(peak_center - tss_pos)
            
#             gene_data[gene_name]['peaks'].add(peak_id)
#             gene_data[gene_name]['total_peak_coverage'] += int(fields[-1])
#             gene_data[gene_name]['distance_to_tss'].append(distance)
        
#         # 6. Create summary DataFrame
#         summary_data = []
#         for gene, data in gene_data.items():
#             if data['peaks']:  # Only include genes with peaks
#                 summary_data.append({
#                     'gene_name': gene,
#                     'num_peaks': len(data['peaks']),
#                     'total_coverage': data['total_peak_coverage'],
#                     'avg_peak_size': data['total_peak_coverage'] / len(data['peaks']),
#                     'min_distance_to_tss': min(data['distance_to_tss']),
#                     'peaks': ','.join(data['peaks'])
#                 })
        
#         df = pd.DataFrame(summary_data)
        
#         # 7. Save results
#         if not df.empty:
#             output_file = os.path.join(output_dir, f"{cell_type}_{condition}_cpg_genes.tsv")
#             df.sort_values('num_peaks', ascending=False).to_csv(output_file, sep='\t', index=False)
#             print(f"Found {len(df)} genes with CpG-overlapping peaks")
#             print(f"Results saved to: {output_file}")
#         else:
#             print(f"No genes found with CpG-overlapping peaks for {cell_type} {condition}")
        
#         return {'genes': df, 'total_peaks': len(cpg_peaks)}
        
#     except Exception as e:
#         print(f"Error processing gene enrichment: {str(e)}")
#         return {'genes': pd.DataFrame(), 'total_peaks': 0}
#     finally:
#         # Cleanup temporary files
#         for f in [tss_bed, extended_tss, temp_peaks]:
#             if os.path.exists(f):
#                 os.remove(f)

# def extract_tss_regions(gtf_file, output_bed):
#     """Extracts TSS regions from GTF file."""
#     with open(gtf_file) as fin, open(output_bed, 'w') as fout:
#         for line in fin:
#             if line.startswith('#'):
#                 continue
            
#             fields = line.strip().split('\t')
#             if fields[2] != "gene":
#                 continue
            
#             # Parse gene information
#             info_dict = {}
#             for item in fields[8].rstrip(';').split('; '):
#                 try:
#                     key, value = item.strip().split(' ', 1)
#                     info_dict[key] = value.strip('"')
#                 except ValueError:
#                     continue
            
#             gene_name = info_dict.get('gene_name', '')
#             gene_type = info_dict.get('gene_type', '')
            
#             if not gene_name or gene_type != "protein_coding":
#                 continue
            
#             # Get TSS position
#             chrom = fields[0]
#             if fields[6] == '+':
#                 tss_pos = int(fields[3])
#                 fout.write(f"{chrom}\t{tss_pos-1}\t{tss_pos}\t{gene_name}\n")
#             else:
#                 tss_pos = int(fields[4])
#                 fout.write(f"{chrom}\t{tss_pos-1}\t{tss_pos}\t{gene_name}\n")

# def create_comparison_summary(results, output_dir):
#     """Creates a summary of gene overlaps between conditions."""
#     summary_file = os.path.join(output_dir, "comparison_summary.txt")
    
#     with open(summary_file, 'w') as f:
#         for cell_type in results:
#             exo_genes = set(results[cell_type]['Exo']['genes']['gene_name']) if not results[cell_type]['Exo']['genes'].empty else set()
#             endo_genes = set(results[cell_type]['Endo']['genes']['gene_name']) if not results[cell_type]['Endo']['genes'].empty else set()
            
#             common_genes = exo_genes & endo_genes
            
#             f.write(f"\n{cell_type} Analysis:\n")
#             f.write(f"Exogenous-specific genes: {len(exo_genes - common_genes)}\n")
#             f.write(f"Endogenous-specific genes: {len(endo_genes - common_genes)}\n")
#             f.write(f"Common genes: {len(common_genes)}\n")
            
#             if common_genes:
#                 common_genes_file = os.path.join(output_dir, f"{cell_type}_common_genes.txt")
#                 with open(common_genes_file, 'w') as cf:
#                     for gene in sorted(common_genes):
#                         cf.write(f"{gene}\n")