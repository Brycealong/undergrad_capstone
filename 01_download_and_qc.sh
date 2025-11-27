#! /usr/bin/bash

#Title: 01 - Data Download and Quality Control
#Purpose: Download raw data from SRA and perform initial FastQC analysis
#To run: nohup bash 01_download_and_qc.sh &

#*************************CONFIGURATION*************************# 

# Project configuration
PROJ_DIR="GSE89692"
BRAIN_REGION="vHIP"
CONDITION="ELS"
OUTPUT_DIR="${PROJ_DIR}/${BRAIN_REGION}/${CONDITION}"

# Directory structure
RAW_DATA="raw_data"
FASTQC_DIR="fastqc"

# Data layout - UPDATE based on your dataset
LAYOUT="paired"      # "paired" or "single"

# Build the SRA list filename dynamically
SRA_LIST_FILE="../../../sra_lists/SRR_Acc_List_${PROJ_DIR}_${BRAIN_REGION}_${CONDITION}.txt"

#*************************SETUP*************************# 

echo "=========================================="
echo "01 - Data Download and Quality Control"
echo "=========================================="
echo "Project: $PROJ_DIR"
echo "Brain Region: $BRAIN_REGION"
echo "Condition: $CONDITION"
echo "Layout: $LAYOUT"
echo "SRA List: $SRA_LIST_FILE"
echo "=========================================="

# Create project directories
echo "Creating output directory: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR || { echo "ERROR: Cannot create directory $OUTPUT_DIR"; exit 1; }
cd $OUTPUT_DIR || { echo "ERROR: Cannot change to directory $OUTPUT_DIR"; exit 1; }
mkdir -p $RAW_DATA $FASTQC_DIR

#*************************DOWNLOAD RAW DATA*************************# 

cd $RAW_DATA

echo ""
echo "[$(date)] Starting SRA data download..."

# Read SRA accessions from file
sra_array=($(cat $SRA_LIST_FILE))
echo "Found ${#sra_array[@]} samples to download"

# Download and convert SRA to FASTQ in one step
# fasterq-dump downloads and converts directly (no need for prefetch)
echo ""
echo "[$(date)] Downloading and converting to FASTQ..."

if [ "$LAYOUT" == "paired" ]; then
	echo "Processing as PAIRED-END data..."
	parallel -j 3 "fasterq-dump {} --split-files --progress" ::: ${sra_array[@]}
elif [ "$LAYOUT" == "single" ]; then
	echo "Processing as SINGLE-END data..."
	parallel -j 3 "fasterq-dump {} --progress" ::: ${sra_array[@]}
else
	echo "ERROR: LAYOUT must be 'paired' or 'single'. Got: $LAYOUT"
	exit 1
fi

# Compress FASTQ files (fasterq-dump outputs uncompressed by default)
echo ""
echo "[$(date)] Compressing FASTQ files..."
if ls *.fastq 1> /dev/null 2>&1; then
	parallel -j 4 gzip {} ::: *.fastq
fi

# Organize files into sample directories
echo ""
echo "[$(date)] Organizing files into sample directories..."
for srr in ${sra_array[@]}; do
	if ls ${srr}*.fastq.gz 1> /dev/null 2>&1; then
		mkdir -p "${srr}"
		mv ${srr}*.fastq.gz ${srr}/
		echo "Organized: $srr"
	else
		echo "WARNING: No FASTQ files found for $srr"
	fi
done

cd ..

#*************************QUALITY CONTROL*************************# 

echo ""
echo "[$(date)] Running FastQC quality control..."

# Run FastQC to quality control and inspect the read length
fastqc -o $FASTQC_DIR -t 10 $RAW_DATA/SRR*/*.fastq.gz

echo ""
echo "=========================================="
echo "[$(date)] Download and QC Complete!"
echo "=========================================="
echo ""
echo "NEXT STEPS:"
echo "1. Review FastQC reports in: $FASTQC_DIR"
echo "2. Check read length from FastQC reports"
echo "3. Verify data quality metrics"
echo "4. Update READ_LENGTH in 02_preprocessing.sh if needed"
echo "5. Ensure LAYOUT='$LAYOUT' in 02_preprocessing.sh matches this run"
echo "6. Run: bash 02_preprocessing.sh"
echo "=========================================="