# MarBTNtranscriptome
## Code to accompany manuscript on gene expression of soft-shell clam transmissible cancer
Samuel Hart, University of Washington Molecular and Cellular Biology PhD program, sfhart33@gmail.com*

<br/>

Follow the steps below to recreate the data analysis from the Metzger lab's pulication on the transcriptomics of transmissible cancer in soft-shell clams. If you use any of these methods or data for your own research, please use the following citation:

* *(In progress... preprint coming fall 2023)*

<br/><br/>

## Dependencies

### Software:
* sratools (v2.10.4) to download RNAseq files
* STAR ()
* blast+ ()

### R (v3.6.0) packages
* tidyverse (v1.3.0)
* msigdbr ()
* DESeq2 ()
* stats ()
* GO.db ()
* factoextra ()
* pheatmap ()
* fgsea ()

## Download github, Mya arenaria genome, and RNAseq files
*Edit target directories and thread count if desired*

Requires [sratools](https://github.com/ncbi/sra-tools/wiki)

Each fastq file may take ~1hr to download - expect to run for a day or two

```
# Download git repo
	git clone https://github.com/sfhart33/MarBTNtranscriptome.git
	cd MarBTNtranscriptome

# Set desired input and ouput folders (output files are minimal, input files are ~400Gb)
	INPUT_FOLDER=/ssd3/RNAseq/inputs
	# INPUT_FOLDER=./inputs
	OUTPUT_FOLDER=/ssd3/RNAseq/outputs
	# OUTPUT_FOLDER=./outputs
	THREADS=50

# Download input data
	./scripts/download_data.sh
```

## Run scripts to prepare STAR index and prepare gene annotations
```
./scripts/index_genome.sh $INPUT_FOLDER $THREADS
./scripts/gene_annotations.sh $INPUT_FOLDER $THREADS
```

## Run makefile to run pipeline
*Edit config.mk with input/output locations and max thread count before running*
```
make

# multithreaded did not work well on our system but may speed things up
	# make --jobs 
	# make --jobs --max-load 50 # 
```
