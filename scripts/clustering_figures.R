#!/usr/bin/Rscript

# Script to plot PCA/clustering figures

# input variables
    output_folder <- commandArgs(trailingOnly=TRUE)[1]
    # output_folder <- "/ssd3/RNAseq/outputs"

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
    names(mycolors) <- sort(unique(annotations$tissue[!(samples_data$tissue %in% c("ASW", "heme_ASW"))]))
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
        rnaseq_norm_tissuespecific[,!(samples_data$tissue %in% c("ASW", "heme_ASW"))] %>% 
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
        include <- !(samples_data$tissue %in% c("ASW", "heme_ASW"))
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
        include <- samples_data$heme %in% c("heme","BTN", "heme_ASW")
        rnaseq_norm_tissuespecific[,include] %>% 
            pheatmap(
            cluster_rows = FALSE,
            show_rownames = FALSE,
            cluster_cols = TRUE,
            clustering_distance_cols = "canberra",
            scale = "row",
            annotation_col = annotations,
            main = "BTN and ASW-BTN and hemocytes, tissue-specific genes"
            )

        rnaseq_norm[rnaseq_norm,include] %>% 
            pheatmap(
            cluster_rows = TRUE,
            show_rownames = FALSE,
            cluster_cols = TRUE,
            clustering_distance_cols = "canberra",
            scale = "row",
            annotation_col = annotations,
            main = "BTN and ASW-BTN and hemocytes, all genes"
            )
dev.off()
pdf("supp_fig4.pdf")
	order <- arrange(rnaseq_norm[(rowSums(rnaseq_norm[,include])> 0),],
				desc(HEL_1_heme)) %>%
				#desc(HEL_1_heme + HEL_2_heme + HEL_3_heme + HEL_4_heme + HEL_5_heme + HEL_6_heme + HEL_7_heme + HEL_8_heme)) %>%
		rownames()
	sub_anno <- structure(list(seq_share = order), .Names = "seq_share", row.names = order, class = "data.frame")
        annotations2 <- annotations
	levels(annotations2$tissue) <- c(levels(annotations2$tissue), "BTN-ASW")
	annotations2$tissue[annotations2$tissue == 'ASW'] <- 'BTN-ASW'
		#annotations2

	rnaseq_norm[order,include] %>% 
            pheatmap(
            cluster_rows = TRUE,
            show_rownames = FALSE,
            cluster_cols = TRUE,
            clustering_distance_cols = "canberra",
            scale = "row",
            annotation_col = annotations2,
            #annotation_row = sub_anno,
            main = "BTN and ASW-BTN and hemocytes, all genes"
            )
dev.off()

# Figure 3

        deseq_vsd <- vst(deseq_asw_btn_heme) # Just ASW and BTN and heme
        fig3 <- plotPCA(deseq_vsd, intgroup = c("tissue")) +
        	theme_classic() +
		theme(legend.position="none") +
        	theme(
            		aspect.ratio = 1,
            		plot.title = element_text(hjust = 0.5),
            		axis.text = element_text(size = 14, face = "bold"),
            		axis.title = element_text(size = 16, face = "bold"),
            		text = element_text(size = 16, face = "bold")
        	)
ggsave("Fig3.pdf", plot = fig3, width = 4, height = 4, units = "in")
ggsave("Fig3.svg", plot = fig3, width = 4, height = 4, units = "in", fix_text_size = FALSE)
