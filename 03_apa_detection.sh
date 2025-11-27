#! /usr/bin/bash

#Title: 03 - APA Site Detection
#Purpose: Detect alternative polyadenylation sites using TAPAS
#To run: nohup bash 03_apa_detection.sh &

#*************************CONFIGURATION*************************# 

# Project configuration
PROJ_DIR="GSE89692"
BRAIN_REGION="vHIP"
CONDITION="ELS"
OUTPUT_DIR="${PROJ_DIR}/${BRAIN_REGION}/${CONDITION}"

# Directory structure
COVERAGE_DIR="coverage"
TAPAS_DIR="tapas"

# Analysis parameters
READ_LENGTH=125  # Must match the read length from your data

#*************************SETUP*************************# 

echo "=========================================="
echo "03 - APA Site Detection"
echo "=========================================="
echo "Project: $PROJ_DIR"
echo "Brain Region: $BRAIN_REGION"
echo "Condition: $CONDITION"
echo "Read Length: $READ_LENGTH"
echo "=========================================="

echo "Using output directory: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR || { echo "ERROR: Cannot create directory $OUTPUT_DIR"; exit 1; }
cd $OUTPUT_DIR || { echo "ERROR: Cannot change to directory $OUTPUT_DIR"; exit 1; }
mkdir -p $TAPAS_DIR

#*************************INSTALL TAPAS*************************# 

echo ""
echo "[$(date)] Setting up TAPAS..."

# Download TAPAS if not already present
if [ ! -d "TAPAS" ]; then
	echo "Cloning TAPAS repository..."
	git clone https://github.com/arefeen/TAPAS.git
fi

# Copy TAPAS tools to working directory
if [ ! -f "APA_sites_detection" ]; then
	echo "Copying TAPAS tools..."
	cp TAPAS/Finding_APA_Sites/* .
	cp TAPAS/Differential_APA_Site_Analysis/* .
	
	# Install required R packages
	if [ -f "package_install.R" ]; then
		echo "Installing R packages..."
		Rscript package_install.R
	fi
fi

echo "[$(date)] TAPAS setup complete."

#*************************APA DETECTION*************************# 

echo ""
echo "[$(date)] Starting APA site detection..."
echo ""

sample_count=0
total_samples=$(ls $COVERAGE_DIR/*.read_coverage.txt | wc -l)

for coverage_file in $COVERAGE_DIR/*.read_coverage.txt; do
	sample_count=$((sample_count + 1))
	
	# Extract sample name
	SAMPLE=$(basename "$coverage_file" .read_coverage.txt)
	
	echo "----------------------------------------"
	echo "[$(date)] Processing sample $sample_count/$total_samples: $SAMPLE"
	echo "----------------------------------------"
	
	# Run APA detection
	./APA_sites_detection \
		-ref mm10.refFlat.txt \
		-cov $coverage_file \
		-l $READ_LENGTH \
		-o $TAPAS_DIR/${SAMPLE}.expression.txt
	
	echo "[$(date)] APA detection complete for $SAMPLE"
	echo ""
done

echo ""
echo "=========================================="
echo "[$(date)] APA Detection Complete!"
echo "=========================================="
echo ""
echo "SUMMARY:"
echo "- Processed $total_samples samples"
echo "- Expression files: $TAPAS_DIR"
echo ""
echo "NEXT STEPS:"
echo "1. Review expression files in $TAPAS_DIR"
echo "2. Configure comparison groups in 04_differential_apa.sh"
echo "3. Run: bash 04_differential_apa.sh"
echo "=========================================="
