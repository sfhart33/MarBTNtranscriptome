# MarBTNtranscriptome
## Code to accompany manuscript on gene expression of soft-shell clam transmissible cancer
Samuel Hart, University of Washington Molecular and Cellular Biology PhD program

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*Contact for questions regarding data analysis: sfhart33@gmail.com*

<br/><br/>

If you use any of these methods or data for your own research, please use the following citation:
*(In progress... preprint coming fall 2023)*

<br/><br/>

## Dependencies

### Software:
* sratools (v2.10.4)
* etc

### R packages
* etc
* etc


## Download Mya arenaria genome and RNAseq files
Edit target directories if desired

Requires [sratools](https://github.com/ncbi/sra-tools/wiki)
```
INPUT_FOLDER=./inputs
OUTPUT_FOLDER=./outputs
mkdir $INPUT_FOLDER/fastq
./scripts/download_data.sh
```

