#!/bin/bash

# Create narrow directory if it doesn't exist
mkdir -p narrow

# Copy and rename each .narrowPeak file
for file in *_peaks.narrowPeak; do
    # Skip if file is a symlink
    if [ -L "$file" ]; then
        continue
    fi
    
    # Extract sample name by removing _peaks.narrowPeak
    sample=${file%_peaks.narrowPeak}
    
    # Copy and rename to narrow directory
    cp "$file" "narrow/${sample}_narrow_peaks.narrowPeak"
done

echo "Files copied and renamed successfully"
