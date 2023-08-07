# Script to plot PCA/clustering figures

# input variables
    output_folder <- commandArgs(trailingOnly=TRUE)[1]

# Load Packages
    library(DESeq2)
    library(tidyverse)
    library(factoextra)
    library(pheatmap)

# Load data tables
    setwd(paste0(output_folder,"/data_tables"))
    load("samples_data.rda")
    load("deseq_untreated.rda")
    load("deseq_asw_pairs.rda")
    load("deseq_asw_alt.rda")
    load("deseq_asw_btn_heme.rda")
    load("deseq_asw_heme.rda")
    load("deseq_run_BTNheme.rda")
    load("deseq_BTNnonheme.rda")
    load("rnaseq.rda")
    load("rnaseq_norm.rda")
    load("rnaseq_norm_tissuespecific.rda")

# MAKE PLOTS
    setwd(paste0(output_folder,"/plots"))

# set color palette and annotations for heir clustering
    cbPalette <- c("#999999",
                   "#E69F00",
                   "#56B4E9",
                   "#009E73",
                   "#F0E442",
                   "#0072B2",
                   "#D55E00",
                   "#CC79A7")
    annotations <- samples_data %>%
        column_to_rownames(var = "new_name") %>%
        select(tissue)
    mycolors <- cbPalette[1:7]
    names(mycolors) <- sort(unique(annotations$tissue[-(4:6)]))
    mycolors <- list(tissue = mycolors)

# Main plots
    pdf("final_pca_and_hc.pdf")
        deseq_vsd <- vst(deseq_untreated) # all non-ASW samples
        plotPCA(deseq_vsd, intgroup = c("tissue")) +
        scale_colour_manual(values=cbPalette) +
        theme_classic() +
        ggtitle("all non-ASW samples") +
        theme(
            aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 16, face = "bold")
        )
        rnaseq_norm_tissuespecific[,-(4:6)] %>% 
            pheatmap(
            cluster_rows = FALSE,
            show_rownames = FALSE,
            cluster_cols = TRUE,
            clustering_distance_cols = "canberra",
            scale = "row",
            annotation_col = annotations,
            annotation_colors = mycolors,
            main = "All samples except ASW, tissue-specific genes"
            )
    dev.off()


# Other plots
pdf("supplementary_pca-hc.pdf")
    # PCA
        deseq_vsd <- vst(deseq_untreated) # all non-ASW samples
        plotPCA(deseq_vsd, intgroup = c("tissue")) +
        theme_classic() +
        ggtitle("all non-ASW samples") +
        theme(
            aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 16, face = "bold")
        )
        deseq_vsd <- vst(deseq_asw_pairs) # Just ASW pairings
        plotPCA(deseq_vsd, intgroup = c("tissue")) +
        theme_classic() +
        ggtitle("Just paired ASW samples") +
        theme(
            aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 16, face = "bold")
        )
        deseq_vsd <- vst(deseq_asw_alt) # Just ASW and BTN
        plotPCA(deseq_vsd, intgroup = c("tissue")) +
        theme_classic() +
        ggtitle("BTN and ASW-BTN") +
        theme(
            aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 16, face = "bold")
        )
        deseq_vsd <- vst(deseq_asw_btn_heme) # Just ASW and BTN and heme
        plotPCA(deseq_vsd, intgroup = c("tissue")) +
        theme_classic() +
        ggtitle("BTN and ASW-BTN and healthy hemocytes") +
        theme(
            aspect.ratio = 1,
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            text = element_text(size = 16, face = "bold")
        )
    # HC   
        rnaseq_norm_tissuespecific[,] %>% 
            pheatmap(
            cluster_rows = FALSE,
            show_rownames = FALSE,
            cluster_cols = TRUE,
            clustering_distance_cols = "canberra",
            scale = "row",
            annotation_col = annotations,
            main = "All samples, tissue-specific genes"
            )
        include <- samples_data$tissue != "ASW"
        rnaseq_norm_tissuespecific[,include] %>% 
            pheatmap(
            cluster_rows = FALSE,
            show_rownames = FALSE,
            cluster_cols = TRUE,
            clustering_distance_cols = "canberra",
            scale = "row",
            annotation_col = annotations,
            main = "All samples except ASW, tissue-specific genes"
            )
        include <- samples_data$heme == "BTN"
        rnaseq_norm_tissuespecific[,include] %>% 
            pheatmap(
            cluster_rows = FALSE,
            show_rownames = FALSE,
            cluster_cols = TRUE,
            clustering_distance_cols = "canberra",
            scale = "row",
            annotation_col = annotations,
            main = "BTN and ASW-BTN, tissue-specific genes"
            )
        rnaseq_norm[,include] %>% 
            pheatmap(
            cluster_rows = FALSE,
            show_rownames = FALSE,
            cluster_cols = TRUE,
            clustering_distance_cols = "canberra",
            scale = "row",
            annotation_col = annotations,
            main = "BTN and ASW-BTN, all genes"
            )
        include <- samples_data$heme %in% c("heme","BTN")
        rnaseq_norm_tissuespecific[,c(1:6,7,13,24,25:28)] %>% 
            pheatmap(
            cluster_rows = FALSE,
            show_rownames = FALSE,
            cluster_cols = TRUE,
            clustering_distance_cols = "canberra",
            scale = "row",
            annotation_col = annotations,
            main = "BTN and ASW-BTN and hemocytes, tissue-specific genes"
            )
        rnaseq_norm[,include] %>% 
            pheatmap(
            cluster_rows = FALSE,
            show_rownames = FALSE,
            cluster_cols = TRUE,
            clustering_distance_cols = "canberra",
            scale = "row",
            annotation_col = annotations,
            main = "BTN and ASW-BTN and hemocytes, all genes"
            )
dev.off()
