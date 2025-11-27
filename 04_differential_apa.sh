#! /usr/bin/bash

#Title: 04 - Differential APA Analysis
#Purpose: Identify differential 3'UTR shortening/lengthening between conditions
#To run: nohup bash 04_differential_apa.sh &

#*************************CONFIGURATION*************************# 

# Project configuration
PROJ_DIR="GSE89692"
BRAIN_REGION="vHIP"
CONDITION="ELS"
OUTPUT_DIR="${PROJ_DIR}/${BRAIN_REGION}/${CONDITION}"

# Directory structure
TAPAS_DIR="tapas"

# Path to metadata file
METADATA_FILE="../../../sample_metadata.csv"

# Define control condition for comparison
CONTROL_CONDITION="Control"
# Treatment condition is automatically set to $CONDITION above

# Analysis parameters
CUTOFF=70           # Cutoff threshold for differential analysis
ANALYSIS_TYPE="s"   

#*************************SETUP*************************# 

echo "=========================================="
echo "04 - Differential APA Analysis"
echo "=========================================="
echo "Project: $PROJ_DIR"
echo "Brain Region: $BRAIN_REGION"
echo "Comparison: $CONTROL_CONDITION vs $CONDITION"
echo "Cutoff: $CUTOFF"
echo "Analysis Type: $ANALYSIS_TYPE"
echo "=========================================="

echo "Using output directory: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR || { echo "ERROR: Cannot create directory $OUTPUT_DIR"; exit 1; }
cd $OUTPUT_DIR || { echo "ERROR: Cannot change to directory $OUTPUT_DIR"; exit 1; }

#*************************BUILD SAMPLE LISTS*************************# 

echo ""
echo "[$(date)] Extracting sample lists from metadata..."

# Extract control samples (C1)
c1_samples=$(awk -F',' -v region="$BRAIN_REGION" -v cond="$CONTROL_CONDITION" \
	'NR>1 && $4==region && $5==cond {print $3}' $METADATA_FILE)

# Extract treatment samples (C2) - using CONDITION variable
c2_samples=$(awk -F',' -v region="$BRAIN_REGION" -v cond="$CONDITION" \
	'NR>1 && $4==region && $5==cond {print $3}' $METADATA_FILE)

# Build comma-separated file paths for C1 (control)
C1_FILES=""
c1_count=0
for sample in $c1_samples; do
	if [ -f "$TAPAS_DIR/${sample}.expression.txt" ]; then
		if [ -z "$C1_FILES" ]; then
			C1_FILES="${TAPAS_DIR}/${sample}.expression.txt"
		else
			C1_FILES="${C1_FILES},${TAPAS_DIR}/${sample}.expression.txt"
		fi
		c1_count=$((c1_count + 1))
	else
		echo "WARNING: Expression file not found for $sample"
	fi
done

# Build comma-separated file paths for C2 (treatment)
C2_FILES=""
c2_count=0
for sample in $c2_samples; do
	if [ -f "$TAPAS_DIR/${sample}.expression.txt" ]; then
		if [ -z "$C2_FILES" ]; then
			C2_FILES="${TAPAS_DIR}/${sample}.expression.txt"
		else
			C2_FILES="${C2_FILES},${TAPAS_DIR}/${sample}.expression.txt"
		fi
		c2_count=$((c2_count + 1))
	else
		echo "WARNING: Expression file not found for $sample"
	fi
done

#*************************VERIFY CONFIGURATION*************************# 

echo ""
echo "=========================================="
echo "Sample Configuration:"
echo "=========================================="
echo "Control ($CONTROL_CONDITION): $c1_count samples"
echo "Treatment ($CONDITION): $c2_count samples"
echo ""
echo "Control files (C1):"
echo "$C1_FILES" | tr ',' '\n' | sed 's/^/  - /'
echo ""
echo "Treatment files (C2):"
echo "$C2_FILES" | tr ',' '\n' | sed 's/^/  - /'
echo "=========================================="

# Check if we have samples for both groups
if [ -z "$C1_FILES" ] || [ -z "$C2_FILES" ]; then
	echo ""
	echo "ERROR: Missing expression files for one or both groups!"
	echo "Please ensure 03_apa_detection.sh completed successfully."
	exit 1
fi

#*************************DIFFERENTIAL ANALYSIS*************************# 

echo ""
echo "[$(date)] Running differential APA analysis..."

# Create output filename based on comparison
OUTPUT_FILE="${CONTROL_CONDITION}_vs_${CONDITION}_${ANALYSIS_TYPE}_result.txt"

# Run differential APA site analysis
./Diff_APA_site_analysis \
	-C1 $C1_FILES \
	-C2 $C2_FILES \
	-a mm10.refFlat.txt \
	-cutoff $CUTOFF \
	-type $ANALYSIS_TYPE \
	-o $OUTPUT_FILE

# Move results to TAPAS directory
if [ -f "diff_result.txt" ]; then
	mv diff_result.txt $TAPAS_DIR/
fi

if [ -f "$OUTPUT_FILE" ]; then
	mv $OUTPUT_FILE $TAPAS_DIR/
fi

echo ""
echo "=========================================="
echo "[$(date)] Differential APA Analysis Complete!"
echo "=========================================="
echo ""
echo "RESULTS:"
echo "- Output directory: $TAPAS_DIR"
echo ""
echo "NEXT STEPS:"
echo "1. Review results in $TAPAS_DIR"
echo "2. Analyze differential APA events"
echo "3. Generate downstream visualizations"
echo "=========================================="
