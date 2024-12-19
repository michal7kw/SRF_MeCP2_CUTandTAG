#!/bin/bash

# Make sure you're in the root of your repository
git filter-repo --force --invert-paths --path 'iterative_alternative/results_5/mecp2_enrichment_independent.csv' \
    --path 'iterative_alternative/results_5/intermediate/endo_combined_peaks.csv' \
    --path 'iterative_alternative/results_5/intermediate/enrichment_results.csv' \
    --path 'iterative_alternative/results_5/intermediate/exo_combined_peaks.csv'

# Force push the changes
git push origin --force 