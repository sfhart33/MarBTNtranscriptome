#!/usr/bin/Rscript

# Script to make volcano plots for deseq and GSEA results

# input variables
    output_folder <- commandArgs(trailingOnly=TRUE)[1]
    # output_folder <- "/ssd3/RNAseq/outputs"

# Load packages
    library(DESeq2)
    library(tidyverse)

# Load deseq and GSEA result data
    setwd(paste0(output_folder,"/data_tables"))
    load("asw_vs_btn_results.rda") # ASW_BTN: negative is higher in ASW
    load("heme_vs_ASWbtn_results.rda") # heme_BTN: negative is higher in BTN
    load("heme_vs_btn_results.rda") # notheme_BTN: negative is higher in BTN
    load("notheme_vs_btn_results.rda") # heme_ASW_BTN: negative is higher in BTN-ASW
    load("ASW_BTN_GSEA.rda")
    load("heme_ASW_BTN_GSEA.rda")
    load("heme_BTN_GSEA.rda")
    load("notheme_BTN_GSEA.rda")

    #new 
    load("UvsP_BTN_results.rda") # UvsP_BTN: negative is higher in PEI
    load("UvsP_heme_results.rda") # UvsP_heme: negative is higher in PEI
    load("heme_vs_aswheme_results.rda") # heme_ASWheme: negative is higher in untreated
    load("USA_vs_PEI_BTN_GSEA.rda")
    load("USA_vs_PEI_heme_GSEA.rda")
    load("heme_vs_ASWheme_GSEA.rda")

# Move to  plots directory
    setwd(paste0(output_folder,"/plots"))

# Volcano plots of all genes for each deseq result
    deseq_volcano <- function(input, comparison){
        sig_limit <- 0.05/nrow(as.data.frame(input))
        # sig_limit <- 0.05
        as.data.frame(input) %>%
            filter(padj < sig_limit) %>%
            ggplot(aes(x = -log2FoldChange, y = -log10(padj))) +
                geom_point(data = filter(as.data.frame(input), padj > sig_limit),
                            aes(x = -log2FoldChange,
                            y = -log10(padj)),
                            color = "grey") +
                geom_point() +
                geom_hline(yintercept =-log10(sig_limit)) +
                xlab(paste("Log2 Fold Change (<-",comparison,"->)")) +
                ylab("-Log10 adj p-value") +
                theme_classic() +
                theme(
                    aspect.ratio = 0.5,
                    plot.title = element_text(hjust = 0.5),
                    axis.text = element_text(size = 14, face = "bold"),
                    axis.title = element_text(size = 16, face = "bold"),
                    text = element_text(size = 16, face = "bold"),
                    legend.title = element_text("none")
                )
    }
    pdf("volcano_plots_genes.pdf")
    deseq_volcano(ASW_BTN, "no treatment BTN vs ASW-BTN")
    deseq_volcano(heme_ASW_BTN, "hemocytes vs ASW-BTN")
    deseq_volcano(heme_BTN, "hemocytes vs BTN")
    deseq_volcano(notheme_BTN, "non-hemocyte tissue vs BTN")

    deseq_volcano(UvsP_BTN, "PEI BTN vs USA BTN")
    deseq_volcano(UvsP_heme, "PEI heme vs USA heme")
    deseq_volcano(heme_ASWheme, "ASW-treated heme vs untreated heme")
    dev.off()

# Volcano plots of all gene gets for each GSEA result
    gsea_volcano <- function(input, comparison){
        sig_limit <- 0.05
        input %>%
            filter(padj < sig_limit) %>%
            ggplot(aes(x = ES, y = -log10(padj))) +
                geom_point(data = filter(input, padj > sig_limit),
                            aes(x = ES,
                            y = -log10(padj)),
                            color = "grey") +
                geom_point() +
                geom_hline(yintercept =-log10(sig_limit)) +
                xlab(paste("Enrichment Score (<-",comparison,"->)")) +
                ylab("-Log10 adj p-value") +
                theme_classic() +
                theme(
                    aspect.ratio = 0.5,
                    plot.title = element_text(hjust = 0.5),
                    axis.text = element_text(size = 14, face = "bold"),
                    axis.title = element_text(size = 16, face = "bold"),
                    text = element_text(size = 16, face = "bold"),
                    legend.title = element_text("none")
                )
    }
    pdf("volcano_plots_gene_sets.pdf")
    gsea_volcano(ASW_BTN_GSEA, "no treatment vs ASW")
    gsea_volcano(heme_ASW_BTN_GSEA, "hemocyte vs ASW-BTN")
    gsea_volcano(heme_BTN_GSEA, "hemocyte vs BTN")
    gsea_volcano(notheme_BTN_GSEA, "non-hemocyte tissue vs BTN")
    gsea_volcano(USA_vs_PEI_BTN_GSEA, "PEI BTN vs USA BTN")
    gsea_volcano(USA_vs_PEI_heme_GSEA, "PEI heme vs USA heme")
    gsea_volcano(heme_vs_ASWheme_GSEA, "ASW-treated heme vs untreated heme")
    dev.off()

# Comparison of results between deseq runs and GSEA runs
    # combine data from deseq runs with directional p-values
    directional_pvalue <- function(input, comparison){
        output <- as.data.frame(input) %>%
            mutate(padj_alt = case_when(
                log2FoldChange > 0 ~ log10(padj),
                log2FoldChange < 0 ~ -log10(padj)
                    )
                ) %>%
            rownames_to_column() %>%
            dplyr::select(rowname, log2FoldChange, padj_alt)
        colnames(output) <- c("gene", paste0(comparison,"_","L2FC"), paste0(comparison,"_","padj"))
        return(output)
    }
    compare_all_genes <- directional_pvalue(ASW_BTN, "BTN_vs_ASW") %>%
        inner_join(directional_pvalue(heme_ASW_BTN, "heme_vs_ASW")) %>%
        inner_join(directional_pvalue(heme_BTN, "heme_vs_BTN")) %>%
        inner_join(directional_pvalue(notheme_BTN, "tissue_vs_BTN")) %>%
        inner_join(directional_pvalue(heme_ASWheme, "ASWheme_vs_heme"))
    # head(compare_all_genes)

    # combine data from GSEA runs with directional p-values
    directional_pvalue_GSEA <- function(input, comparison){
        output <- input %>%
            mutate(padj_alt = case_when(
                NES > 0 ~ -log10(padj),
                NES < 0 ~ log10(padj)
                    )
                ) %>%
            dplyr::select(pathway, NES, padj_alt)
        colnames(output) <- c("pathway", paste0(comparison,"_","NES"), paste0(comparison,"_","padj"))
        return(output)
    }
    compare_all_gene_sets <- directional_pvalue_GSEA(ASW_BTN_GSEA, "BTN_vs_ASW") %>%
        inner_join(directional_pvalue_GSEA(heme_ASW_BTN_GSEA, "heme_vs_ASW")) %>%
        inner_join(directional_pvalue_GSEA(heme_BTN_GSEA, "heme_vs_BTN")) %>%
        inner_join(directional_pvalue_GSEA(notheme_BTN_GSEA, "tissue_vs_BTN")) %>%
        inner_join(directional_pvalue_GSEA(heme_vs_ASWheme_GSEA,  "ASWheme_vs_heme"))
    # head(compare_all_gene_sets)

    # Function to make comparison plots
    directional_pvalue_plots <- function(x_axis, y_axis, x_label, y_label, type){
        if(!type %in% c("deseq", "gsea") ) {stop("Must be deseq or gsea input")}
        if(type == "deseq") {
            sig_limit <- 0.05/43039
            data_set <- compare_all_genes
        }
        if(type == "gsea") {
            sig_limit <- 0.05
            data_set <- compare_all_gene_sets
        }
        ggplot(data_set,
               aes(get(paste0(x_axis,"_","padj")),
                   get(paste0(y_axis,"_","padj")))) +
            geom_point() +
            geom_hline(yintercept = -log10(sig_limit)) +
            geom_hline(yintercept = log10(sig_limit)) +
            geom_vline(xintercept = -log10(sig_limit)) +
            geom_vline(xintercept = log10(sig_limit)) +
            xlab(paste(x_label, "\nadjusted p-value (Log10, +/- directional)")) +
            ylab(paste(y_label, "\nadjusted p-value (Log10, +/- directional)")) +
            theme_classic() +
            theme(
                aspect.ratio = 1,
                plot.title = element_text(hjust = 0.5),
                axis.text = element_text(size = 14, face = "bold"),
                axis.title = element_text(size = 16, face = "bold"),
                text = element_text(size = 16, face = "bold"),
                legend.title = element_text("none")
            )
    }
    # deseq gene plots
    pdf("comparison_plots_genes.pdf")
    directional_pvalue_plots("heme_vs_BTN", "tissue_vs_BTN", "BTN vs hemocytes", "BTN vs tissue", "deseq")
    directional_pvalue_plots("heme_vs_BTN", "heme_vs_ASW", "BTN vs hemocytes", "ASW-BTN vs hemocytes", "deseq")
    directional_pvalue_plots("BTN_vs_ASW", "heme_vs_ASW", "ASW-BTN vs BTN", "ASW-BTN vs hemocytes", "deseq")
    directional_pvalue_plots("BTN_vs_ASW", "heme_vs_BTN", "ASW-BTN vs BTN", "BTN vs hemocytes", "deseq")
    directional_pvalue_plots("BTN_vs_ASW", "ASWheme_vs_heme", "ASW-BTN vs BTN", "hemocytes vs ASW-heme", "deseq")
    dev.off()
    # GSEA gene set plots
    pdf("comparison_plots_gene_sets.pdf")
    directional_pvalue_plots("heme_vs_BTN", "tissue_vs_BTN", "BTN vs hemocytes", "BTN vs tissue", "gsea")
    directional_pvalue_plots("heme_vs_BTN", "heme_vs_ASW", "BTN vs hemocytes", "ASW-BTN vs hemocytes", "gsea")
    directional_pvalue_plots("BTN_vs_ASW", "heme_vs_ASW", "ASW-BTN vs BTN", "ASW-BTN vs hemocytes", "gsea")
    directional_pvalue_plots("BTN_vs_ASW", "heme_vs_BTN", "ASW-BTN vs BTN", "BTN vs hemocytes", "gsea")
    directional_pvalue_plots("BTN_vs_ASW", "ASWheme_vs_heme", "ASW-BTN vs BTN", "hemocytes vs ASW-heme", "gsea")
    dev.off()