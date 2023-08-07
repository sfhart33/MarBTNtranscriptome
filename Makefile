# Run pipeline for M arenaria bivalve transmissible neoplasia RNAseq analysis

# Variable definitions
include config.mk

FASTQS=$(wildcard $(INPUT_FOLDER)/fastq/*R1_001.fastq.gz)
STAROUT=$(patsubst $(INPUT_FOLDER)/fastq/%_R1_001.fastq.gz, $(OUTPUT_FOLDER)/star/%_ReadsPerGene.out.tab, $(FASTQS))
DESEQOUT=$(awk -v location="$(OUTPUT_FOLDER)/star/" '{print location $0}' deseq_output.lst)
# First target that checks if output files are made
.PHONY : outputs
outputs : $(OUTFILES)

# Test print
#$(OUTPUT_FOLDER)/test/%_ReadsPerGene.out.tab : $(INPUT_FOLDER)/fastq/%_R1_001.fastq.gz $(INPUT_FOLDER)/fastq/%_R2_001.fastq.gz
#	echo $^ $* > $@


# 1) Run STAR alignment
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

# 2) Run DEseq on reads/gene data
.PHONY : nextstep
$(DESEQOUT) : $(STAROUT)
	mkdir -p $(OUTPUT_FOLDER)/star