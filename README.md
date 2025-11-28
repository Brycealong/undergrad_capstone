# RNA-seq APA Analysis Pipeline

A modular pipeline for analyzing alternative polyadenylation (APA) in RNA-seq data using TAPAS.

## Pipeline Overview

This pipeline is split into 4 independent, runnable scripts:

1. **01_download_and_qc.sh** - Download data from SRA and perform quality control
2. **02_preprocessing.sh** - Trim adapters, align to genome, generate coverage files
3. **03_apa_detection.sh** - Detect APA sites using TAPAS
4. **04_differential_apa.sh** - Identify differential APA between conditions

## Quick Start

### 1. Configure Your Analysis

Edit the configuration section at the top of each script:

```bash
PROJ_DIR="GSE89692"
BRAIN_REGION="vHIP"        # vHIP, mPFC, NAc, or VTA
CONDITION="ELS"             # ELS, CSDS, or EC
```

For script 04, also configure your comparison:

```bash
CONTROL_CONDITION="Control"
# CONDITION variable (set above) is automatically used as treatment
```

Note: The treatment condition is already defined by the `CONDITION` variable in scripts 01-03. Script 04 compares Control vs whatever condition you processed.

### 2. Run the Pipeline

```bash
# Step 1: Download and QC
nohup bash 01_download_and_qc.sh > 01_download.log 2>&1 &

# After completion, review FastQC reports
# Check: GSE89692/vHIP/ELS/fastqc/

# Step 2: Preprocessing and alignment
# Update READ_LENGTH in script 02 based on FastQC results
nohup bash 02_preprocessing.sh > 02_preprocessing.log 2>&1 &

# Step 3: APA detection
nohup bash 03_apa_detection.sh > 03_apa_detection.log 2>&1 &

# Step 4: Differential APA analysis
nohup bash 04_differential_apa.sh > 04_differential_apa.log 2>&1 &
```

## Directory Structure

```
GSE89692/
└── vHIP/
    └── ELS/
        ├── raw_data/          # Raw FASTQ files
        ├── fastqc/            # Quality control reports
        ├── trimmed/           # Adapter-trimmed reads
        ├── sam/               # Alignment files (SAM)
        ├── bam/               # Sorted BAM files
        ├── coverage/          # Read coverage files
        └── tapas/             # APA expression and results
```

## Required Files

### Before Running:

1. **SRA accession list**: `sra_lists/SRR_Acc_List_GSE89692_vHIP_ELS.txt`
   - One SRR accession per line
   - Automatically constructed from: `SRR_Acc_List_${PROJ_DIR}_${BRAIN_REGION}_${CONDITION}.txt`
   - you can obtain other SRR runs at: https://www.ncbi.nlm.nih.gov/Traces/study/
2. **Metadata file**: `sample_metadata.csv`
   - Used for grouping samples in differential analysis
   - metadata: https://docs.google.com/spreadsheets/d/1tusEfvgCp4mUn4E8Y3s2XV0iERSjC7R5PibH4COlb3o/edit?usp=sharing. You should download it as a csv file and rename the file to `sample_metadata.csv`

## Important Parameters

### Script 01: Download and QC
- `parallel -j 3`: Number of parallel downloads (adjust based on connection)

### Script 02: Preprocessing
- `READ_LENGTH`: **MUST UPDATE** after reviewing FastQC results from script 01
- `NUM_THREADS`: Adjust based on available CPU cores

### Script 03: APA Detection
- `READ_LENGTH`: Must match your actual read length

### Script 04: Differential Analysis
- `CUTOFF`: Threshold for differential analysis (default: 70)

## Checkpoints Between Scripts

### After Script 01:
✓ Review FastQC HTML reports in `fastqc/`
✓ Note the read length from FastQC
✓ Check for adapter contamination
✓ Verify data quality scores

### After Script 02:
✓ Check alignment rates from HISAT2 logs
✓ Verify BAM files exist and are sorted
✓ Check coverage file sizes

### After Script 03:
✓ Verify expression files exist for all samples
✓ Check file sizes are reasonable (not empty)

## Troubleshooting

### Script 01 fails to download
- Check internet connection
- Verify SRA accessions in list file
- Ensure SRA toolkit is installed

### Script 02 fails during alignment
- Verify `READ_LENGTH` is set correctly
- Check that reference genome downloaded completely
- Ensure sufficient disk space

### Script 03 produces empty expression files
- Verify coverage files are not empty
- Check that `READ_LENGTH` matches your data
- Ensure mm10.refFlat.txt exists

### Script 04 reports missing files
- Ensure script 03 completed successfully
- Verify sample IDs in metadata match expression files
- Check that metadata CSV format is correct

## Dependencies

- [ncbi/sra-tools](https://github.com/ncbi/sra-tools) (fasterq-dump)
- [samtools](https://www.htslib.org/)
- [R](https://www.r-project.org/) (for TAPAS)
- [GNU parallel](https://www.gnu.org/software/parallel/) 
- [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
- [HISAT2](https://daehwankimlab.github.io/hisat2/)
- FastQC
- Git
- Standard Unix tools (awk, grep, sed)

## Output Files

### Main Results (in `tapas/` directory):

- `{SAMPLE}.expression.txt` - APA expression for each sample
- `Control_vs_ELS_s_result.txt` - Differential shortening APA results
- `diff_result.txt` - Detailed differential analysis

## Tips

1. **Disk Space**: RNA-seq analysis requires significant storage. Monitor disk usage.
2. **Parallel Jobs**: Adjust `-j` parameter based on your system resources.
3. **Memory**: Large BAM files may require substantial RAM for sorting.
5. **Rerunning Scripts**: Each script can be run independently. If a step fails, fix the issue and rerun just that script.

## Citation

If using this pipeline, please cite:
- TAPAS: Ashraful Arefeen, Juntao Liu, Xinshu Xiao, Tao Jiang, TAPAS: tool for alternative polyadenylation site analysis, *Bioinformatics*, Volume 34, Issue 15, August 2018, Pages 2521–2529, https://doi.org/10.1093/bioinformatics/bty110
- HISAT2: Kim, D., Paggi, J.M., Park, C. *et al.* Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nat Biotechnol* **37**, 907–915 (2019). https://doi.org/10.1038/s41587-019-0201-4
- other appropriate tools if necessary

## Support

For issues or questions, check:
- TAPAS documentation: https://github.com/arefeen/TAPAS
- FastQC manual
- HISAT2 manual

  

  
