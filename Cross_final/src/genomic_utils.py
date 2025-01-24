"""Utility functions for genomic operations."""

import pandas as pd
import numpy as np
from typing import List, Tuple, Dict
import pybedtools
import warnings

def load_gtf(gtf_file: str) -> pd.DataFrame:
    """Load gene coordinates from GTF file."""
    print(f"Loading gene annotations from {gtf_file}")
    genes = []
    
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'gene':
                # Parse attributes
                attrs = {}
                for attr in fields[8].split(';'):
                    if not attr.strip():
                        continue
                    try:
                        key, value = attr.strip().split(' ', 1)
                        attrs[key] = value.strip('"')
                    except ValueError:
                        continue
                
                # Get gene name, trying different possible attribute names
                gene_id = (attrs.get('gene_name') or 
                          attrs.get('gene_id') or 
                          attrs.get('Name'))
                
                if gene_id:
                    genes.append({
                        'gene': gene_id,
                        'chromosome': fields[0],
                        'start': int(fields[3]),
                        'end': int(fields[4]),
                        'strand': fields[6],
                        'tss': int(fields[3]) if fields[6] == '+' else int(fields[4])
                    })
                    
    print(f"Loaded {len(genes)} genes from GTF")
    return pd.DataFrame(genes).set_index('gene')

def get_promoter_coords(gene_coords: pd.DataFrame, 
                       window: Tuple[int, int]) -> pd.DataFrame:
    """Get promoter coordinates for genes."""
    promoters = gene_coords.copy()
    promoters['start'] = promoters.apply(
        lambda x: x['tss'] + window[0] if x['strand'] == '+' 
        else x['tss'] - window[1], axis=1
    )
    promoters['end'] = promoters.apply(
        lambda x: x['tss'] + window[1] if x['strand'] == '+' 
        else x['tss'] - window[0], axis=1
    )
    return promoters

def get_gene_body_coords(gene_coords: pd.DataFrame, 
                        extension: int) -> pd.DataFrame:
    """Get gene body coordinates with extension."""
    gene_bodies = gene_coords.copy()
    gene_bodies['start'] = gene_bodies.apply(
        lambda x: x['start'] - extension if x['strand'] == '+' 
        else x['start'], axis=1
    )
    gene_bodies['end'] = gene_bodies.apply(
        lambda x: x['end'] if x['strand'] == '+' 
        else x['end'] + extension, axis=1
    )
    return gene_bodies

def calculate_peak_overlaps(peaks: List[str], 
                          regions: pd.DataFrame,
                          threshold: float = 0.25) -> pd.DataFrame:
    """Calculate overlap between peaks and genomic regions."""
    # Convert regions to BedTool format
    regions_df = regions.reset_index()
    regions_df = regions_df[['chromosome', 'start', 'end', 'gene']]
    regions_bed = pybedtools.BedTool.from_dataframe(regions_df)
    
    overlaps = []
    for peak_file in peaks:
        try:
            peaks_bed = pybedtools.BedTool(peak_file)
            # Use -wo to get both original entries (A and B) plus the overlap width
            intersect = regions_bed.intersect(peaks_bed, wo=True)
            
            for hit in intersect:
                region_length = hit.end - hit.start
                # The overlap width is the last field when using -wo
                overlap_length = float(hit[-1])
                if overlap_length / region_length >= threshold:
                    overlaps.append({
                        'gene': hit.name,
                        'peak_file': peak_file,
                        'overlap_ratio': overlap_length / region_length
                    })
        except Exception as e:
            print(f"Error processing peak file {peak_file}: {str(e)}")
            continue
    
    return pd.DataFrame(overlaps)

def normalize_bigwig_values(values: List[float], 
                          method: str = 'rpm') -> List[float]:
    """Normalize bigwig values."""
    if method == 'rpm':
        total = sum(v for v in values if v is not None)
        if total > 0:
            return [v * 1e6 / total if v is not None else 0 for v in values]
    return [v if v is not None else 0 for v in values]
