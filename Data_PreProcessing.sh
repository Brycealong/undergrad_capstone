#! usr/bin/sh

#Title: Data downloading and pre-processing
#To run - nohup sh Data_PreProcessing.sh

#*************************Downloading Raw Data[SRA]*************************# 

## Download raw data from SRA using GNU parallel and SRA-toolkit

#Create project directory to store raw data
PROJ_DIR="GSE89692"
BRAIN_REGION="vHIP"
CONTROL_CONDITION="Control"
TREATMENT_CONDITION="ELS"
READ_LEN=125
OUTPUT_DIR="${PROJ_DIR}/${BRAIN_REGION}/${TREATMENT_CONDITION}"
RAW_DATA="raw_data"
FASTQC_DIR="fastqc"
TRIM_DIR="trimmed"
SAM_DIR="sam"
BAM_DIR="bam"
COVERAGE_DIR="coverage"
TAPAS_DIR="tapas"

# Path to metadata file
METADATA_FILE="../../../sample_metadata.csv"

# Build the SRA list filename dynamically
SRA_LIST_FILE="../../../sra_lists/SRR_Acc_List_${PROJ_DIR}_${BRAIN_REGION}_${TREATMENT_CONDITION}.txt"

cd $OUTPUT_DIR
mkdir -p $RAW_DATA $FASTQC_DIR $TRIM_DIR $SAM_DIR $BAM_DIR $COVERAGE_DIR $TAPAS_DIR
cd $RAW_DATA

echo "Using SRA list file: $SRA_LIST_FILE"

# Download SRA files using accession numbers
parallel -j 3 prefetch {} ::: $(cat $SRA_LIST_FILE)

#Convert SRA to fastq
parallel -j 3 fasterq-dump {} ::: $(cat $SRA_LIST_FILE)

# Get a list of unique SRR numbers
srr_numbers=$(ls | grep -o 'SRR[0-9]\+' | sort -u)

# Create directories for each SRR number
for srr_number in ${srr_numbers}; do
	mkdir -p "${srr_number}"

	# Move files to the corresponding directory
	mv ${srr_number}*.fastq.gz ${srr_number}/
done

cd ..

# fastqc to quality control and inspect the read length
fastqc -o $FASTQC_DIR -t 10 $RAW_DATA/SRR*/*.fastq.gz # threads 10

#*************************PRE-PROCESSING*************************# 
# # Unzip the fastq.gz data to gain fastq files
# parallel "gunzip {}" ::: $RAW_DATA/*.fastq.gz

# download the adapter file
curl -O https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa

# Download the HISAT2 index of mm10
curl -O https://genome-idx.s3.amazonaws.com/hisat/mm10_genome.tar.gz
tar -vxzf mm10_genome.tar.gz
# the indexes will be extracted to mm10/

# Download the refFlat format of mm10
curl -O https://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/refFlat.txt.gz
unzip refFlat.txt.gz
mv refFlat.txt mm10.refFlat.txt #rename

# Download the TAPAS tool
git clone https://github.com/arefeen/TAPAS.git
cp Finding_APA_Sites/* .
cp Differential_APA_Site_Analysis/* .
Rscript package_install.R script

for dir in $RAW_DATA/SRR*; do
	# Find the R1 and R2 FASTQ files
	r1_file=$(find "$dir" -name "*_1.fastq.gz")
	r2_file=$(find "$dir" -name "*_2.fastq.gz")

	# Extract the sample name
	OUTPUT=$(basename "$dir")

	# Apply adapter trimming
	trimmomatic PE -threads 12 "$r1_file" "$r2_file" \
    	"$TRIM_DIR/${OUTPUT}_1.trimmed.fastq.gz" \
    	"$TRIM_DIR/${OUTPUT}_1un.trimmed.fastq.gz" \
    	"$TRIM_DIR/${OUTPUT}_2.trimmed.fastq.gz" \
    	"$TRIM_DIR/${OUTPUT}_2un.trimmed.fastq.gz" \
    	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

	# Perform genome alignment of reads using HISAT2
	hisat2 -p 12 -x mm10/genome \
		-1 $TRIM_DIR/${OUTPUT}_1.trimmed.fastq.gz \
		-2 $TRIM_DIR/${OUTPUT}_2.trimmed.fastq.gz \
		-S $SAM_DIR/"${OUTPUT}.sam"

	# Convert sam files to bam files
	samtools view -b $SAM_DIR/${OUTPUT}.sam -o $BAM_DIR/${OUTPUT}.bam
	samtools sort $BAM_DIR/${OUTPUT}.bam -o $BAM_DIR/${OUTPUT}.sorted.bam
	samtools index $BAM_DIR/${OUTPUT}.sorted.bam

	# Convert bam files to read coverage files
	samtools depth $BAM_DIR/"${OUTPUT}.sorted.bam" -o $COVERAGE_DIR/"${OUTPUT}.read_coverage.txt"	

	# Perform APA detection
	./APA_sites_detection -ref mm10_refFlat.txt \
		-cov $COVERAGE_DIR/${OUTPUT}.read_coverage.txt \
		-l 125 -o $TAPAS_DIR/${OUTPUT}.expression.txt
done

# Perform shortening/lengthening detection
./Diff_APA_site_analysis \
	-C1 $TAPAS_DIR/SRR8956213.expression.txt,$TAPAS_DIR/SRR8956225.expression.txt,$TAPAS_DIR/SRR8956253.expression.txt,$TAPAS_DIR/SRR8956257.expression.txt,$TAPAS_DIR/SRR8956265.expression.txt,$TAPAS_DIR/SRR8956293.expression.txt \
	-C2 $TAPAS_DIR/SRR8956214.expression.txt,$TAPAS_DIR/SRR8956226.expression.txt,$TAPAS_DIR/SRR8956254.expression.txt,$TAPAS_DIR/SRR8956258.expression.txt,$TAPAS_DIR/SRR8956266.expression.txt,$TAPAS_DIR/SRR8956294.expression.txt \
	-a mm10_refFlat.txt -cutoff 70 -type s \
	-o shortening_result_final_tapas.txt

mv diff_result.txt shortening_result_final_tapas.txt $TAPAS_DIR

#*************************DIFFERENTIAL APA ANALYSIS*************************# 

# Extract samples from metadata based on brain region and conditions
# Using awk to parse CSV and filter by brain_region and stress_condition
c1_samples=$(awk -F',' -v region="$BRAIN_REGION" -v cond="$CONTROL_CONDITION" \
	'NR>1 && $4==region && $5==cond {print $3}' $METADATA_FILE)

c2_samples=$(awk -F',' -v region="$BRAIN_REGION" -v cond="$TREATMENT_CONDITION" \
	'NR>1 && $4==region && $5==cond {print $3}' $METADATA_FILE)

# Build comma-separated lists with full paths to expression files
C1_FILES=""
for sample in $c1_samples; do
	if [ -z "$C1_FILES" ]; then
		C1_FILES="${TAPAS_DIR}/${sample}.expression.txt"
	else
		C1_FILES="${C1_FILES},${TAPAS_DIR}/${sample}.expression.txt"
	fi
done

C2_FILES=""
for sample in $c2_samples; do
	if [ -z "$C2_FILES" ]; then
		C2_FILES="${TAPAS_DIR}/${sample}.expression.txt"
	else
		C2_FILES="${C2_FILES},${TAPAS_DIR}/${sample}.expression.txt"
	fi
done

# Perform shortening/lengthening detection
./Diff_APA_site_analysis \
	-C1 $C1_FILES \
	-C2 $C2_FILES \
	-a mm10_refFlat.txt -cutoff 70 -type s \
	-o shortening_result_final_tapas.txt

mv diff_result.txt shortening_result_final_tapas.txt $TAPAS_DIR

