#!/usr/bin/Rscript

# Script to run DEseq comparisons and plot PCA/clustering figures

# input variables
    input_folder <- commandArgs(trailingOnly=TRUE)[1]
    # input_folder <- "/ssd3/RNAseq/inputs"
	# output_folder <- "/ssd3/RNAseq/outputs"
    output_folder <- commandArgs(trailingOnly=TRUE)[2]

# Load Packages
    library(DESeq2)
    library(tidyverse)

# Load sample information file 
#NOTE NAMES MAY BE DIFFERENT IF DOWNLOADED FROM SRA
    print("Loading sample data")
    # samples_data <- read.table("/ssd3/RNAseq/STAR_genome_output/28samples/sample_data.txt", header = TRUE) %>%
    samples_data <- read.table("./inputs/sample_data.txt", header = TRUE) %>%
        mutate(new_name = paste(clam, tissue, sep = "_"))
    samples <- samples_data %>% pull(name) %>% as.vector()

# Healthy subset
    samples_data_healthy <- samples_data %>%
        filter(clam %in% c("HEL_1", "HEL_2", "HEL_3"))
    samples_healthy <- samples_data_healthy %>% pull(name) %>% as.vector()

# Gene name conversion (these are in same order)
    #gene_names <- read.table("/ssd1/mya_genome_ncbi_annotation/gene_translation.txt",
    gene_names <- paste0(input_folder, "/GCF_026914265.1_ASM2691426v1/gene_translation.txt") %>%
        read.table(sep="\t", header = FALSE, col.names = c("mya_gene","desc"), quote = "")
    head(gene_names)

# Load data   
    # setwd("/ssd3/RNAseq/STAR_NCBI_genome_output")
    setwd(paste0(output_folder, "/star"))
    for(sample in samples){
        file <- paste0(sample, "_ReadsPerGene.out.tab")
        file_load <- read.table(file, sep="\t", skip = 4,
                                col.names = c("mya_gene",
                                sample, "fwd", "rev")) %>%
            select(-fwd,-rev)
        if(sample == samples[1]){
            print("first")
            rnaseq <- file_load 
        } else {
            print(sample)
            rnaseq <- full_join(rnaseq,file_load)
        }
    }
    # rename with better columns
        rnaseq <- column_to_rownames(rnaseq, var = "mya_gene")
        colnames(rnaseq) <- samples_data$new_name
        head(rnaseq)
        rnaseq_healthy <- rnaseq[,as.vector(samples_data_healthy$new_name)]

    # Normalize rnaseq data frame
        norm_value <- max(colSums(rnaseq))/colSums(rnaseq)
        rnaseq_norm <- as.data.frame(t(t(as.matrix(rnaseq))*norm_value))


###################
# DESEQ COMPARISONS
###################
    print("Running DEseq comparisons")
    #setwd("/home/shart/github/MarBTNtranscriptome/outputs/data_tables")
    setwd(paste0(output_folder, "/data_tables"))
    # Deseq on everything except treated cells
        include <- !(samples_data$tissue %in% c("ASW", "heme_ASW"))
        deseq_untreated <- DESeqDataSetFromMatrix(countData = rnaseq[,include],
                                    colData = samples_data[include,],
                                    design = ~ tissue) %>%
            DESeq()
        save(deseq_untreated, file = "deseq_untreated.rda")

    # Deseq on just ASW pairings
        include <- samples_data$clam %in% c("BTN_1", "BTN_2", "BTN_3")
        deseq_asw_pairs <- DESeqDataSetFromMatrix(countData = rnaseq[,include],
                                    colData = samples_data[include,],
                                    design = ~ clam + tissue) %>%
            DESeq()
        save(deseq_asw_pairs, file = "deseq_asw_pairs.rda")
        ASW_BTN <- results(deseq_asw_pairs, name="tissue_BTN_vs_ASW")
        save(ASW_BTN, file = "asw_vs_btn_results.rda")

    # Deseq on all BTN vs ASW
        include <- (samples_data$heme == "BTN")
        deseq_asw_alt <- DESeqDataSetFromMatrix(countData = rnaseq[,include],
                                    colData = samples_data[include,],
                                    design = ~ tissue) %>%
            DESeq()
        save(deseq_asw_alt, file = "deseq_asw_alt.rda")

    # Deseq on ASW and BTN and heme
        include <- samples_data$heme %in% c("BTN", "heme", "heme_ASW")
        deseq_asw_btn_heme <- DESeqDataSetFromMatrix(countData = rnaseq[,include],
                                    colData = samples_data[include,],
                                    design = ~ tissue) %>%
            DESeq()
        save(deseq_asw_btn_heme, file = "deseq_asw_btn_heme.rda")

    # Deseq on just ASW vs heme
        include <- samples_data$tissue %in% c("ASW", "heme")
        deseq_asw_heme <- DESeqDataSetFromMatrix(countData = rnaseq[,include],
                                    colData = samples_data[include,],
                                    design = ~ tissue) %>%
            DESeq()
        save(deseq_asw_heme, file = "deseq_asw_heme.rda")
        heme_ASW_BTN <- results(deseq_asw_heme, name="tissue_heme_vs_ASW")
        save(heme_ASW_BTN, file = "heme_vs_ASWbtn_results.rda")

    # Deseq on BTN vs hemocytes
        include <- samples_data$tissue %in% c("BTN", "heme")
        deseq_run_BTNheme <- DESeqDataSetFromMatrix(countData = rnaseq[,include],
                                    colData = samples_data[include,],
                                    design = ~ tissue) %>%
            DESeq()
        save(deseq_run_BTNheme, file = "deseq_run_BTNheme.rda")
        heme_BTN <- results(deseq_run_BTNheme, name="tissue_heme_vs_BTN")
        save(heme_BTN, file = "heme_vs_btn_results.rda")

    # Deseq on BTN vs everything but hemocytes
        include <- !(samples_data$tissue %in% c("ASW", "heme", "heme_ASW"))
        deseq_BTNnonheme <- DESeqDataSetFromMatrix(countData = rnaseq[,include],
                                    colData = samples_data[include,],
                                    design = ~ heme) %>%
            DESeq()
        save(deseq_BTNnonheme, file = "deseq_BTNnonheme.rda")
        notheme_BTN <- results(deseq_BTNnonheme, name="heme_not_vs_BTN")
        save(notheme_BTN, file = "notheme_vs_btn_results.rda")

            #############################
            #### ADD PEI COMPARISONS ####
            #############################
#    # Deseq on BTN USA vs PEI
#        include <- (samples_data$tissue == "BTN")
#        deseq_UvsP_BTN <- DESeqDataSetFromMatrix(countData = rnaseq[,include],
#                                    colData = samples_data[include,],
#                                    design = ~ location) %>%
#            DESeq()
#        save(deseq_UvsP_BTN, file = "deseq_UvsP_BTN.rda")
#        UvsP_BTN <- results(deseq_UvsP_BTN, name="location_USA_vs_PEI")
#        save(UvsP_BTN, file = "UvsP_BTN_results.rda")

    # Deseq on healthy USA vs PEI
        include <- (samples_data$tissue == "heme")
        deseq_UvsP_heme <- DESeqDataSetFromMatrix(countData = rnaseq[,include],
                                    colData = samples_data[include,],
                                    design = ~ location) %>%
            DESeq()
        save(deseq_UvsP_heme, file = "deseq_UvsP_heme.rda")
        UvsP_heme <- results(deseq_UvsP_heme, name="location_USA_vs_PEI")
        save(UvsP_heme, file = "UvsP_heme_results.rda")

    # Deseq on heme ASW pairings
        include <- samples_data$clam %in% c("HEL_6", "HEL_7", "HEL_8")
        deseq_asw_heme_pairs <- DESeqDataSetFromMatrix(countData = rnaseq[,include],
                                    colData = samples_data[include,],
                                    design = ~ clam + tissue) %>%
            DESeq()
        save(deseq_asw_heme_pairs, file = "deseq_asw_heme_pairs.rda")
        heme_ASWheme <- results(deseq_asw_heme_pairs, name="tissue_heme_ASW_vs_heme")
        save(heme_ASWheme, file = "heme_vs_aswheme_results.rda")
###############################################
# DESEQ TO DETERMINE TISSUE-SPECIFIC EXPRESSION
###############################################
    print("Tissue-by-tissue DEseq comparisons to determine tissue-specific genes")
    # Blank df
        tissue_specific_top100 <- data.frame(gene = c(), tissue = c())
    #tissue list
        tissues <- samples_data_healthy$tissue %>%
            unique() %>%
            as.vector()

    for(i in tissues){
    # new column specififying target
        samples_data_healthy2 <- samples_data_healthy %>%
            mutate(target = get(i))
    # run deseq and extract top 100 genes expressed in each tissue
        tissue_specific_top100 <- DESeqDataSetFromMatrix(
                                    countData = rnaseq_healthy,
                                    colData = samples_data_healthy2,
                                    design = ~ clam + target) %>%
            DESeq() %>%
            results(contrast=c("target","not",i)) %>%
            as.data.frame() %>%
            arrange(stat) %>%
            mutate(tissue = i) %>%
            rownames_to_column(var = "gene") %>%
            select(gene, tissue) %>%
            head(n=100) %>%
            rbind(tissue_specific_top100)
    }
    tissue_specific_list <- unique(tissue_specific_top100$gene)

# normalized values for tissue-specific list
    rnaseq_norm_tissuespecific <- rnaseq_norm[tissue_specific_list,]

    write.table(rnaseq, file = "rnaseq.tsv", quote = FALSE, sep = "\t")
    save(samples_data, file = "samples_data.rda")
    save(rnaseq, file = "rnaseq.rda")
    save(rnaseq_norm, file = "rnaseq_norm.rda")
    save(rnaseq_norm_tissuespecific, file = "rnaseq_norm_tissuespecific.rda")
    print("Done with DEseq")