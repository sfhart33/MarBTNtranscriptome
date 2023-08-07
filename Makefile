# Run pipeline for M arenaria bivalve transmissible neoplasia RNAseq analysis

# Variable definitions
include config.mk

FASTQS=$(wildcard $(INPUT_FOLDER)/fastq/*R1_001.fastq.gz)
STAROUT=$(patsubst $(INPUT_FOLDER)/fastq/%_R1_001.fastq.gz, $(OUTPUT_FOLDER)/star/%_ReadsPerGene.out.tab, $(FASTQS))
DESEQOUT=$(shell awk -v location=$(OUTPUT_FOLDER)/data_tables/ '{print location $$0}' deseq_output.lst)
GSEAOUT=$(shell awk -v location=$(OUTPUT_FOLDER)/data_tables/ '{print location $$0}' gsea_output.lst)
CLUSTFIGS=$(shell awk -v location=$(OUTPUT_FOLDER)/plots/ '{print location $$0}' clustering_figs.lst)
VOLCANOFIGS=$(shell awk -v location=$(OUTPUT_FOLDER)/plots/ '{print location $$0}' volcano_figs.lst)

# First target that checks if output files are made
.PHONY : outputs
outputs : $(CLUSTFIGS) $(VOLCANOFIGS)

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
$(DESEQOUT) : $(STAROUT)
	Rscript ./scripts/deseq_comparisons.R $(INPUT_FOLDER) $(OUTPUT_FOLDER)

# 3) Run GSEA
$(GSEAOUT) : $(DESEQOUT)
	Rscript ./scripts/GSEA.R $(INPUT_FOLDER) $(OUTPUT_FOLDER)

# 4) Print sample clustering figures
$(CLUSTFIGS) : $(DESEQOUT)
	Rscript ./scripts/clustering_figures.R $(OUTPUT_FOLDER)

# 5) Print volcano plot figures
$(VOLCANOFIGS) : $(DESEQOUT) $(GSEAOUT)
	Rscript ./scripts/volcano_plots.R $(OUTPUT_FOLDER)