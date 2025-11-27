# undergrad_capstone

## Dependencies

- [ncbi/sra-tools](https://github.com/ncbi/sra-tools)
- [samtools](https://www.htslib.org/)
- [R](https://www.r-project.org/)
- [GNU parallel](https://www.gnu.org/software/parallel/) 
- [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)
- [HISAT2](https://daehwankimlab.github.io/hisat2/)

```
conda install -c bioconda sra-tools samtools trimmomatic hisat2
conda install -c conda-forge r-essentials r-base parallel
```

The spreadsheet of the metadata of this study is in: https://docs.google.com/spreadsheets/d/1tusEfvgCp4mUn4E8Y3s2XV0iERSjC7R5PibH4COlb3o/edit?usp=sharing

### Step 1: Data Preprocessing

```
nohup sh Data_PreProcessing.sh
```

