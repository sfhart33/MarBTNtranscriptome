# MarBTNtranscriptome
## Code to accompany manuscript on gene expression of soft-shell clam transmissible cancer
Samuel Hart, University of Washington Molecular and Cellular Biology PhD program

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;*Contact for questions regarding data analysis: sfhart33@gmail.com*

Follow the steps below to recreate the data analysis from the Metzger lab's pulication on the transcriptomics of transmissible cancer in soft-shell clams. If you use any of these methods or data for your own research, please use the following citation:

* *(In progress... preprint coming fall 2023)*

<br/><br/>

## Dependencies

### Software:
* sratools (v2.10.4) to download RNAseq files
* STAR
* blast+

### R packages
* etc
* etc


## Download github, Mya arenaria genome, and RNAseq files
*Edit target directories and thread count if desired*

Requires [sratools](https://github.com/ncbi/sra-tools/wiki)
```
git clone https://github.com/sfhart33/MarBTNtranscriptome.git
cd MarBTNtranscriptome

INPUT_FOLDER=/ssd3/RNAseq/inputs
# INPUT_FOLDER=./inputs
OUTPUT_FOLDER=/ssd3/RNAseq/outputs
# OUTPUT_FOLDER=./outputs
THREADS=50

mkdir $INPUT_FOLDER/fastq

./scripts/download_data.sh
```

## Run scripts to prepare STAR index and prepare gene annotations
```
./scripts/index_genome.sh $INPUT_FOLDER $THREADS
./scripts/gene_annotations.sh $INPUT_FOLDER $THREADS
```

## Run makefile to run pipeline
```
make # standard
# make --jobs # multithreaded
#make --jobs --max-load 50 # multithreaded with max threads
```
