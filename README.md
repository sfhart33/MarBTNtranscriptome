# MarBTNtranscriptome
## Code to accompany manuscript on gene expression of soft-shell clam transmissible cancer
Samuel Hart, University of Washington Molecular and Cellular Biology PhD program, sfhart33@gmail.com*

<br/>

Follow the steps below to recreate the data analysis from the Metzger lab's pulication on the transcriptomics of transmissible cancer in soft-shell clams. If you use any of these methods or data for your own research, please use the following citation:

* *(In progress... preprint coming fall 2023)*

<br/><br/>

## Dependencies

### Software:
* sratools (2.10.4) to download RNAseq files
* STAR (2.7.4a) to align reads
* blast (2.10.0) to annotate closest gene hits
* singularity (3.6.4) to run...
  * STAR-Fusion (1.11.0) to call gene fusions
  * Can alter code to run other STAR-Fusion installations

### R (3.6.0) packages used
* tidyverse (1.3.0)
* msigdbr (7.5.1)
* DESeq2 (1.26.0)
* GO.db (3.10.0)
* factoextra (1.0.7)
* pheatmap (1.0.12)
* fgsea (1.27.0)

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

# Create file structure
	./scripts/mkdirs.sh

# Download input data
	./scripts/download_data.sh
```

## Run scripts to prepare STAR index and prepare gene annotations
```
./scripts/index_genome.sh $INPUT_FOLDER $THREADS
./scripts/gene_annotations.sh $INPUT_FOLDER $THREADS
```

## Run makefile to run main pipeline
*Edit config.mk with input/output locations and max thread count before running*
```
make outputs
```

## For fusion analysis 
*We used singularity image of STAR-Fusion. You may need to alter these scripts depending on how you access STAR-Fusion.*
```
# Download singularity image
	wget https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/__archived_releases/star-fusion.v1.11.0.simg -P $INPUT_FOLDER/

# Generate genome index to run STAR-Fusion
	./scripts/index_starfusion.sh $INPUT_FOLDER $THREADS

# Run Fusion portion of pipeline
	make fusion
```

