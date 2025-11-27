#! /usr/bin/bash

#Title: 02 - Data Preprocessing and Alignment
#Purpose: Trim adapters, align to genome, and generate coverage files
#To run: nohup bash 02_preprocessing.sh &

#*************************CONFIGURATION*************************# 

# Project configuration
PROJ_DIR="GSE89692"
BRAIN_REGION="vHIP"
CONDITION="ELS"
OUTPUT_DIR="${PROJ_DIR}/${BRAIN_REGION}/${CONDITION}"

# Directory structure
RAW_DATA="raw_data"
TRIM_DIR="trimmed"
SAM_DIR="sam"
BAM_DIR="bam"
COVERAGE_DIR="coverage"

# Analysis parameters
LAYOUT="paired"      # "paired" or "single" - UPDATE based on your data
READ_LENGTH=125      # UPDATE THIS based on FastQC results from script 01
NUM_THREADS=12       # Adjust based on your system

#*************************SETUP*************************# 

echo "=========================================="
echo "02 - Data Preprocessing and Alignment"
echo "=========================================="
echo "Project: $PROJ_DIR"
echo "Brain Region: $BRAIN_REGION"
echo "Condition: $CONDITION"
echo "Layout: $LAYOUT"
echo "Read Length: $READ_LENGTH"
echo "Threads: $NUM_THREADS"
echo "=========================================="

echo "Using output directory: $OUTPUT_DIR"
mkdir -p $OUTPUT_DIR || { echo "ERROR: Cannot create directory $OUTPUT_DIR"; exit 1; }
cd $OUTPUT_DIR || { echo "ERROR: Cannot change to directory $OUTPUT_DIR"; exit 1; }
mkdir -p $TRIM_DIR $SAM_DIR $BAM_DIR $COVERAGE_DIR

#*************************DOWNLOAD REFERENCE FILES*************************# 

echo ""
echo "[$(date)] Downloading reference files and tools..."

# Download adapter files if not present
if [ "$LAYOUT" == "paired" ] && [ ! -f "TruSeq3-PE.fa" ]; then
	echo "Downloading Trimmomatic paired-end adapters..."
	curl -O https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa
elif [ "$LAYOUT" == "single" ] && [ ! -f "TruSeq3-SE.fa" ]; then
	echo "Downloading Trimmomatic single-end adapters..."
	curl -O https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-SE.fa
fi

# Download HISAT2 index if not present
if [ ! -d "mm10" ]; then
	echo "Downloading HISAT2 mm10 index..."
	curl -O https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
	tar -vxzf mm10_genome.tar.gz
	echo "HISAT2 index extracted to mm10/"
fi

# Download refFlat annotation if not present
if [ ! -f "mm10.refFlat.txt" ]; then
	echo "Downloading mm10 refFlat annotation..."
	curl -O https://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refFlat.txt.gz
	gunzip refFlat.txt.gz
	mv refFlat.txt mm10.refFlat.txt
fi

echo "[$(date)] Reference files ready."

#*************************PREPROCESSING PIPELINE*************************# 

echo ""
echo "[$(date)] Starting preprocessing pipeline..."
echo ""

sample_count=0
total_samples=$(ls -d $RAW_DATA/SRR* | wc -l)

for dir in $RAW_DATA/SRR*; do
	sample_count=$((sample_count + 1))
	
	# Extract the sample name
	SAMPLE=$(basename "$dir")
	
	echo "----------------------------------------"
	echo "[$(date)] Processing sample $sample_count/$total_samples: $SAMPLE"
	echo "----------------------------------------"

	if [ "$LAYOUT" == "paired" ]; then
		#*************************PAIRED-END PROCESSING*************************#
		
		# Find the R1 and R2 FASTQ files
		r1_file=$(find "$dir" -name "*_1.fastq.gz")
		r2_file=$(find "$dir" -name "*_2.fastq.gz")
		
		if [ -z "$r1_file" ] || [ -z "$r2_file" ]; then
			echo "ERROR: Could not find paired-end files for $SAMPLE"
			echo "Expected: ${dir}/*_1.fastq.gz and ${dir}/*_2.fastq.gz"
			continue
		fi
		
		echo "[$(date)] Found R1: $(basename $r1_file)"
		echo "[$(date)] Found R2: $(basename $r2_file)"

		# Step 1: Adapter trimming (Paired-End)
		echo "[$(date)] Trimming adapters (PE mode)..."
		trimmomatic PE -threads $NUM_THREADS "$r1_file" "$r2_file" \
			"$TRIM_DIR/${SAMPLE}_1.trimmed.fastq.gz" \
			"$TRIM_DIR/${SAMPLE}_1un.trimmed.fastq.gz" \
			"$TRIM_DIR/${SAMPLE}_2.trimmed.fastq.gz" \
			"$TRIM_DIR/${SAMPLE}_2un.trimmed.fastq.gz" \
			ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		# Step 2: Genome alignment (Paired-End)
		echo "[$(date)] Aligning to genome (PE mode)..."
		hisat2 -p $NUM_THREADS -x mm10/genome \
			-1 $TRIM_DIR/${SAMPLE}_1.trimmed.fastq.gz \
			-2 $TRIM_DIR/${SAMPLE}_2.trimmed.fastq.gz \
			-S $SAM_DIR/"${SAMPLE}.sam"
	
	elif [ "$LAYOUT" == "single" ]; then
		#*************************SINGLE-END PROCESSING*************************#
		
		# Find the single FASTQ file
		r1_file=$(find "$dir" -name "*.fastq.gz" | head -n 1)
		
		if [ -z "$r1_file" ]; then
			echo "ERROR: Could not find FASTQ file for $SAMPLE"
			echo "Expected: ${dir}/*.fastq.gz"
			continue
		fi
		
		echo "[$(date)] Found: $(basename $r1_file)"

		# Step 1: Adapter trimming (Single-End)
		echo "[$(date)] Trimming adapters (SE mode)..."
		trimmomatic SE -threads $NUM_THREADS "$r1_file" \
			"$TRIM_DIR/${SAMPLE}.trimmed.fastq.gz" \
			ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

		# Step 2: Genome alignment (Single-End)
		echo "[$(date)] Aligning to genome (SE mode)..."
		hisat2 -p $NUM_THREADS -x mm10/genome \
			-U $TRIM_DIR/${SAMPLE}.trimmed.fastq.gz \
			-S $SAM_DIR/"${SAMPLE}.sam"
	
	else
		echo "ERROR: LAYOUT must be 'paired' or 'single'. Got: $LAYOUT"
		exit 1
	fi

	# Step 3: SAM to BAM conversion and sorting (same for both)
	echo "[$(date)] Converting to BAM and sorting..."
	samtools view -b $SAM_DIR/${SAMPLE}.sam -o $BAM_DIR/${SAMPLE}.bam
	samtools sort $BAM_DIR/${SAMPLE}.bam -o $BAM_DIR/${SAMPLE}.sorted.bam
	samtools index $BAM_DIR/${SAMPLE}.sorted.bam
	
	# Optional: Remove intermediate files to save space
	# rm $SAM_DIR/${SAMPLE}.sam
	# rm $BAM_DIR/${SAMPLE}.bam

	# Step 4: Generate coverage file (same for both)
	echo "[$(date)] Generating read coverage..."
	samtools depth $BAM_DIR/"${SAMPLE}.sorted.bam" -o $COVERAGE_DIR/"${SAMPLE}.read_coverage.txt"
	
	echo "[$(date)] Sample $SAMPLE complete!"
	echo ""
done

echo ""
echo "=========================================="
echo "[$(date)] Preprocessing Complete!"
echo "=========================================="
echo ""
echo "SUMMARY:"
echo "- Processed $total_samples samples"
echo "- Trimmed reads: $TRIM_DIR"
echo "- Aligned BAM files: $BAM_DIR"
echo "- Coverage files: $COVERAGE_DIR"
echo ""
echo "NEXT STEPS:"
echo "1. Review alignment statistics"
echo "2. Check BAM file quality"
echo "3. Run: bash 03_apa_detection.sh"
echo "=========================================="
