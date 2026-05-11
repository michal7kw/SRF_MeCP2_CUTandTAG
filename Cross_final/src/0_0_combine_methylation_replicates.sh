#!/bin/bash

# --- Configuration ---
# Path to the chromosome sizes file
CHROM_SIZES="../../data/mm10.chrom.sizes"
# Directory containing the input BigWig files
# DATA_DIR=".../data/bigwig"
DATA_DIR="../results_local/normalized_bigwig"

# Directory to save the merged BigWig files
OUTPUT_DIR="../../data/bigwig/merged_normalized" # Changed to a subdirectory

# --- Tool Check ---
command -v wiggletools >/dev/null 2>&1 || { echo >&2 "Error: wiggletools not found. Please install it (e.g., via conda or https://github.com/Ensembl/WiggleTools)."; exit 1; }
command -v wigToBigWig >/dev/null 2>&1 || { echo >&2 "Error: wigToBigWig not found. Please install UCSC Kent tools (http://hgdownload.soe.ucsc.edu/admin/exe/)."; exit 1; }

# --- Function to merge replicates ---
merge_bw_files() {
  local output_prefix="$1"
  shift # Remove output_prefix from arguments
  local input_files=("$@") # Remaining arguments are input files

  local temp_wig="${OUTPUT_DIR}/${output_prefix}.wig"
  local final_bw="${OUTPUT_DIR}/${output_prefix}.bw"
  local inputs_str=""

  # Check if all input files exist
  for file in "${input_files[@]}"; do
    if [ ! -f "$file" ]; then
      echo "Error: Input file not found: $file"
      return 1 # Indicate error
    fi
    inputs_str+=" $file" # Build string for wiggletools command
  done

  echo "Merging to ${final_bw}..."

  # Calculate mean using wiggletools, output to temp Wiggle file
  wiggletools mean ${inputs_str} > "$temp_wig"

  # Convert Wiggle to BigWig using wigToBigWig
  wigToBigWig "$temp_wig" "$CHROM_SIZES" "$final_bw"

  # Clean up temporary Wiggle file
  rm "$temp_wig"

  echo "Created ${final_bw}"
  return 0 # Indicate success
}

# --- Main Execution ---
echo "Starting BigWig replicate merging..."

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# # Merge DP Input replicates (r1, r2, r3)
# merge_bw_files \
#   "Medip_DP_input_merged" \
#   "${DATA_DIR}/Medip_DP_input_r1.bw" \
#   "${DATA_DIR}/Medip_DP_input_r2.bw" \
#   "${DATA_DIR}/Medip_DP_input_r3.bw"

# # Merge DP Output replicates (r1, r3)
# merge_bw_files \
#   "Medip_DP_output_merged" \
#   "${DATA_DIR}/Medip_DP_output_r1.bw" \
#   "${DATA_DIR}/Medip_DP_output_r3.bw"

# # Merge N Input replicates (r1, r2, r3)
# merge_bw_files \
#   "Medip_N_input_merged" \
#   "${DATA_DIR}/Medip_N_input_r1.bw" \
#   "${DATA_DIR}/Medip_N_input_r2.bw" \
#   "${DATA_DIR}/Medip_N_input_r3.bw"

# # Merge N Output replicates (r2, r3)
# merge_bw_files \
#   "Medip_N_output_merged" \
#   "${DATA_DIR}/Medip_N_output_r2.bw" \
#   "${DATA_DIR}/Medip_N_output_r3.bw"

# # Merge PP Input replicates (r1, r2, r3)
# merge_bw_files \
#   "Medip_PP_input_merged" \
#   "${DATA_DIR}/Medip_PP_input_r1.bw" \
#   "${DATA_DIR}/Medip_PP_input_r2.bw" \
#   "${DATA_DIR}/Medip_PP_input_r3.bw"

# # Merge PP Output replicates (r1, r3) - No change needed here
# merge_bw_files \
#   "Medip_PP_output_merged" \
#   "${DATA_DIR}/Medip_PP_output_r1.bw" \
#   "${DATA_DIR}/Medip_PP_output_r3.bw"


########################

# Merge DP Output replicates (r1, r3)
merge_bw_files \
  "Medip_DP_output_merged" \
  "${DATA_DIR}/DP_r1_methylation_signal.bw" \
  "${DATA_DIR}/DP_r3_methylation_signal.bw"

# Merge N Output replicates (r2, r3)
merge_bw_files \
  "Medip_N_output_merged" \
  "${DATA_DIR}/N_r2_methylation_signal.bw" \
  "${DATA_DIR}/N_r3_methylation_signal.bw"

# Merge PP Output replicates (r1, r3) - No change needed here
merge_bw_files \
  "Medip_PP_output_merged" \
  "${DATA_DIR}/PP_r1_methylation_signal.bw" \
  "${DATA_DIR}/PP_r3_methylation_signal.bw"

echo "-----------------------------"
echo "Merging process completed."
echo "Merged files are in: ${OUTPUT_DIR}"
echo "-----------------------------"