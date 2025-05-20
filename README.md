# undergrad_capstone

## Dependencies

- [ncbi/sra-tools](https://github.com/ncbi/sra-tools)
- [samtools](https://www.htslib.org/)
- [R](https://www.r-project.org/)
- [GNU parallel](https://www.gnu.org/software/parallel/) 
- [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)

```
conda install -c bioconda sra-tools samtools trimmomatic
conda install -c conda-forge r-essentials r-base parallel
```

### Step 1: Data Preprocessing

```
nohup sh Data_PreProcessing.sh
```

