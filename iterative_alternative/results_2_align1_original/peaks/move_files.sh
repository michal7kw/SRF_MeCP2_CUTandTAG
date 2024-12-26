#!/bin/bash

# Create narrow directory if it doesn't exist
mkdir -p narrow

# Move all narrowPeak files to the narrow directory
cp *_narrow_peaks.narrowPeak narrow/

echo "Copied all narrowPeak files to ./narrow directory"
