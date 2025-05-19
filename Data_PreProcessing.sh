#! usr/bin/sh

#Title: Data downloading and pre-processing
#To run - nohup sh Data_PreProcessing.sh

#*************************Downloading Raw Data[SRA]*************************# 

## Download raw data from SRA using GNU parallel and SRA-toolkit

#Creat project directory to store raw data
PROJ_DIR="GSE89692"
BRAIN_REGION="vHIP"
CONDITION="ELS"
OUTPUT_DIR="${PROJ_DIR}/${BRAIN_REGION}/${CONDITION}"
RAW_DATA="raw_data"
FASTQC_DIR="fastqc"
TRIM_DIR="trimmed"
ALIGNED_DIR="aligned"
BAM_DIR="bam"
COVERAGE_DIR="coverage"
TAPAS_DIR="tapas"

cd $OUTPUT_DIR
mkdir -p $RAW_DATA $FASTQC_DIR $TRIM_DIR $ALIGNED_DIR $BAM_DIR $COVERAGE_DIR $TAPAS_DIR
cd $RAW_DATA
# Download SRA files using accession numbers
parallel -j 3 prefetch {} ::: \
SRR8956213 SRR8956225 SRR8956253 SRR8956257 SRR8956265 SRR8956293 \
SRR8956214 SRR8956226 SRR8956254 SRR8956258 SRR8956266 SRR8956294

#Convert SRA to fastq
parallel -j 3 fastq-dump --gzip --skip-technical --read-filter pass --dumpbase --split-3 --clip --origfmt {} ::: \
SRR8956213/SRR8956213.sra SRR8956225/SRR8956225.sra SRR8956253/SRR8956253.sra SRR8956257/SRR8956257.sra SRR8956265/SRR8956265.sra SRR8956293/SRR8956293.sra \
SRR8956214/SRR8956214.sra SRR8956226/SRR8956226.sra SRR8956254/SRR8956254.sra SRR8956258/SRR8956258.sra SRR8956266/SRR8956266.sra SRR8956294/SRR8956294.sra

cd ..

#*************************PRE-PROCESSING*************************# 
# # Unzip the fastq.gz data to gain fastq files
# parallel "gunzip {}" ::: $RAW_DATA/*.fastq.gz

# download the adapter file
curl -O https://raw.githubusercontent.com/usadellab/Trimmomatic/main/adapters/TruSeq3-PE.fa

# Apply adapter trimming
for dir in $RAW_DATA/SRR*; do
	# Find the R1 and R2 FASTQ files
	r1_file=$(find "$dir" -name "*_1.fastq.gz")
	r2_file=$(find "$dir" -name "*_2.fastq.gz")

	# Extract the sample name
	samp=$(basename "$dir")

	trimmomatic PE -threads 12 "$r1_file" "$r2_file" \
    	"${TRIM_DIR}/${samp}_1.trimmed.fastq.gz" \
    	"${TRIM_DIR}/${samp}_1un.trimmed.fastq.gz" \
    	"${TRIM_DIR}/${samp}_2.trimmed.fastq.gz" \
    	"${TRIM_DIR}/${samp}_2un.trimmed.fastq.gz" \
    	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done
