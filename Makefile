# Run pipeline for M arenaria bivalve transmissible neoplasia RNAseq analysis

# Variable definitions
include config.mk

FASTQS=$(wildcard $(INPUT_FOLDER)/fastq/*R1_001.fastq.gz)
#FASTQS=$(INPUT_FOLDER)/fastq/21-FFM-15G1_R1_001.fastq.gz
OUTFILES=$(patsubst $(INPUT_FOLDER)/fastq/%_R1_001.fastq.gz, $(OUTPUT_FOLDER)/star/%_ReadsPerGene.out.tab, $(FASTQS))

# First target that checks if output files are made
.PHONY : outputs
outputs : $(OUTFILES)

# Test print
#$(OUTPUT_FOLDER)/test/%_ReadsPerGene.out.tab : $(INPUT_FOLDER)/fastq/%_R1_001.fastq.gz $(INPUT_FOLDER)/fastq/%_R2_001.fastq.gz
#	echo $^ $* > $@


# Run STAR alignment
$(OUTPUT_FOLDER)/star/%_ReadsPerGene.out.tab : $(INPUT_FOLDER)/fastq/%_R1_001.fastq.gz $(INPUT_FOLDER)/fastq/%_R2_001.fastq.gz
	STAR \
            --runThreadN $(THREADS) \
            --genomeDir $(INPUT_FOLDER)/STARindex \
            --readFilesIn $^ \
            --outFileNamePrefix $(OUTPUT_FOLDER)/star/$*"_" \
            --readFilesCommand zcat \
            --quantMode GeneCounts \
            --outSAMtype None
	rm $(OUTPUT_FOLDER)/star/$*_SJ.out.tab