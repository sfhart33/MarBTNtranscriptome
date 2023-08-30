# Run pipeline for M arenaria bivalve transmissible neoplasia RNAseq analysis

# Variable definitions
include config.mk

FASTQS=$(wildcard $(INPUT_FOLDER)/fastq/*R1_001.fastq.gz)
STAROUT=$(patsubst $(INPUT_FOLDER)/fastq/%_R1_001.fastq.gz, $(OUTPUT_FOLDER)/star/%_ReadsPerGene.out.tab, $(FASTQS))
FUSIONOUT=$(patsubst $(INPUT_FOLDER)/fastq/%_R1_001.fastq.gz, $(OUTPUT_FOLDER)/star_fusion/%.fusion_predictions.abridged.tsv, $(FASTQS))
DESEQOUT=$(shell awk -v location=$(OUTPUT_FOLDER)/data_tables/ '{print location $$0}' deseq_output.lst)
GSEAOUT=$(shell awk -v location=$(OUTPUT_FOLDER)/data_tables/ '{print location $$0}' gsea_output.lst)
CLUSTFIGS=$(shell awk -v location=$(OUTPUT_FOLDER)/plots/ '{print location $$0}' clustering_figs.lst)
VOLCANOFIGS=$(shell awk -v location=$(OUTPUT_FOLDER)/plots/ '{print location $$0}' volcano_figs.lst)

# Command to make output files for main pipeline
.PHONY : outputs
outputs : $(CLUSTFIGS) $(VOLCANOFIGS) # $(GSEAOUT)

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

# Command to run fusion analysis
.PHONY : fusion
fusion : $(FUSIONOUT)

# 6) Run STAR-Fusion and save output file
$(OUTPUT_FOLDER)/star_fusion/%.fusion_predictions.abridged.tsv : $(INPUT_FOLDER)/fastq/%_R1_001.fastq.gz $(INPUT_FOLDER)/fastq/%_R2_001.fastq.gz
	echo $*
	singularity exec -e \
		-B `pwd`,/ssd3/RNAseq/fusions_ncbi/ctat_genome_lib_build_dir,$(INPUT_FOLDER)/fastq,$(OUTPUT_FOLDER)/star_fusion \
		/ssd3/RNAseq/fusions/star-fusion.v1.11.0.simg \
		STAR-Fusion \
			--left_fq  $(INPUT_FOLDER)/fastq/$*"_R1_001.fastq.gz" \
			--right_fq $(INPUT_FOLDER)/fastq/$*"_R2_001.fastq.gz" \
			--genome_lib_dir /ssd3/RNAseq/fusions_ncbi/ctat_genome_lib_build_dir \
			--output_dir $(OUTPUT_FOLDER)/star_fusion/$*"_StarFusionOut" \
			--FusionInspector validate \
			--examine_coding_effect \
			--denovo_reconstruct \
			--CPU $(THREADS)
	cp $(OUTPUT_FOLDER)/star_fusion/$*"_StarFusionOut/star-fusion.fusion_predictions.abridged.tsv" $(OUTPUT_FOLDER)/star_fusion/$*".fusion_predictions.abridged.tsv"
	rm -r $(OUTPUT_FOLDER)/star_fusion/$*"_StarFusionOut"

# Command to run genome effects analysis
.PHONY : genome_effects
genome_effects : $(OUTPUT_FOLDER)/plots/fusion_plot.pdf $(OUTPUT_FOLDER)/plots/cn_exp_plot.pdf

# 7) Analyze fusion gene overlap among samples and how copy number and Steamer insertions affect expression
$(OUTPUT_FOLDER)/plots/fusion_plot.pdf $(OUTPUT_FOLDER)/plots/cn_exp_plot.pdf : $(FUSIONOUT) ./inputs/genes_by_PEI_copy_number.tsv ./inputs/genes_by_USA_copy_number.tsv ./inputs/steamer_genes.tsv
	Rscript ./scripts/gene_effects_analysis.R $(INPUT_FOLDER) $(OUTPUT_FOLDER)