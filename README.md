# MarBTNtranscriptome
## Code to accompany manuscript on gene expression of soft-shell clam (*Mya arenaria*) bivalve transmissible neoplasia (MarBTN) 

Samuel Hart, University of Washington Molecular and Cellular Biology PhD program, *sfhart33@gmail.com*

<br/>

Follow the steps below to recreate the data analysis from the [Metzger lab](https://pnri.org/metzger-lab/)'s publication on the transcriptomics of transmissible cancer in soft-shell clams. If you use any of these methods or data for your own research, please use the following citation:

* Manuscript currently under review... (preprint will be posted Sept 2024)
  * *Gene expression in soft-shell clam (*Mya arenaria*) transmissible cancer reveals survival mechanisms during host infection and seawater transfer* 
* In the meantime, preliminary results and discussion can be found in the second chapter of my dissertation (results as of May 2023)
  * [*Evolution of soft-shell clam transmissible cancer*](https://digital.lib.washington.edu/researchworks/handle/1773/50503)

<br/><br/>

## Dependencies

### Software:
* [sratools](https://github.com/ncbi/sra-tools/wiki) (2.10.4) to download RNAseq files
* [STAR](https://github.com/alexdobin/STAR) (2.7.4a) to align reads
* [blast+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#downloadblastdata) (2.10.0) to annotate closest gene hits
* [singularity](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html) (3.6.4) to run...
  * [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion) (1.11.0) to call gene fusions
  * Can alter code to run other STAR-Fusion installations

### R (3.6.0) packages:
* tidyverse (1.3.0) to code and plot
* msigdbr (7.5.1) to annotate gene sets
* GO.db (3.10.0) to annotate gene sets
* factoextra (1.0.7) for PCA
* pheatmap (1.0.12) for clustered heatmaps
* DESeq2 (1.26.0) for differential expression analysis
* fgsea (1.27.0) for gene set enrichment analysis

## Download github, Mya arenaria genome, and RNAseq files
*Edit target directories and thread count if desired*

Requires [sratools](https://github.com/ncbi/sra-tools/wiki)

Each fastq file may take ~1hr to download - expect to run for a day or two

*To run downstream scripts, use NCBI SRA metadata files to rename SRA fastq files (e.g. SRR23856981_1.fastq.gz) to the original name (e.g. 01-FFM-28E5_R1_001.fastq.gz). Or email sfhart33@gmail.com: I plan to incorperate a translation script but have not yet done so.*

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

Warning: STAR-Fusion can take a LONG time to run for all samples
```
# Download singularity image
	wget https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/STAR-Fusion/__archived_releases/star-fusion.v1.11.0.simg -P $INPUT_FOLDER/

# Generate genome index to run STAR-Fusion
	./scripts/index_starfusion.sh $INPUT_FOLDER $THREADS

# Run Fusion portion of pipeline
	make fusion
```

Alternatively, simply use outputs we have already included on github to save time by running the following
```
cp -i ./inputs/star_fusion $OUTPUT_FOLDER/star_fusion
```

## To analyze the effect of unstable genome on gene expression
```
make genome_effects
```

## To confirm BTN samples are from USA sub-lineage 

notes in ./scripts/mt_snvs.txt

* Use STAR to align RNA reads to mt genome
* Checked manually in IGV for SNPs
* Compare to USA/PEI SNPS FROM GENOME PAPER
  * External: https://github.com/sfhart33/MarBTNgenome/tree/main/06_Mito_analysis
  * Internal Metzger lab: /ssd3/Mar_genome_analysis/bwa_mapping/mito/all_samples/somatypus/Somatypus_SNVs_final.counts
* Use SNP ratios as rough proxy for purity too
