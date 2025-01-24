"""Configuration settings for the cross-analysis pipeline."""

import os

# Base directories
BASE_DIR = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_CUTandTAG"
MEDIP_DIR = "/beegfs/scratch/ric.broccoli/kubacki.michal/DATA/MECP2/MEDIP"
SMARCB1_DIR = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1/results"
RESULTS_DIR = os.path.join(BASE_DIR, "Cross_final/results")
DATA_DIR = os.path.join(BASE_DIR, "Cross_final/data")

# Reference files
GTF_FILE = os.path.join(BASE_DIR, "DATA/gencode.vM10.annotation.gtf")

# Input files
GENE_LISTS = {
    'up': os.path.join(DATA_DIR, 'up.csv'),
    'down': os.path.join(DATA_DIR, 'down.csv'),
    'no_deg': os.path.join(DATA_DIR, 'no_deg.csv')
}

# MeDIP data paths
MEDIP_BIGWIG = {
    'N': [
        os.path.join(MEDIP_DIR, f'output_done/bigwig/Medip_N_output_r{i}.bw')
        for i in range(1, 4)
    ],
    'PP': [
        os.path.join(MEDIP_DIR, f'output_done/bigwig/Medip_PP_output_r{i}.bw')
        for i in range(1, 4)
    ],
    'DP': [
        os.path.join(MEDIP_DIR, f'output_done/bigwig/Medip_DP_output_r{i}.bw')
        for i in range(1, 4)
    ]
}

# SMARCB1 CUT&Tag data
SMARCB1_PEAKS = {
    'BG': [
        os.path.join(SMARCB1_DIR, f'peaks_alt/BG{i}_peaks.narrowPeak')
        for i in range(1, 4)
    ],
    'BM': [
        os.path.join(SMARCB1_DIR, 'peaks_alt/BM3_peaks.narrowPeak')
    ]
}

# Analysis parameters
PROMOTER_WINDOW = (-2000, 2000)  # Window around TSS
GENE_BODY_EXTENSION = 1000  # Extend gene body by this many bp
PEAK_OVERLAP_THRESHOLD = 0.25  # Minimum overlap required for peak intersection
