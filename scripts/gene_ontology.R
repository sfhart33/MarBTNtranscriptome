#!/usr/bin/Rscript

# Script to merge genes with ontologies for GSEA downstream

# Load packages
        library(tidyverse)
        library(msigdbr)


# Load data
        setwd(commandArgs(trailingOnly=TRUE))
        uniprot <- read.delim("GCF_026914265.1_ASM2691426v1_closest_gene_name.txt", sep="\t", head=FALSE, col.names = c("mRNA_gene","uniprot_hit"))
        mRNAs <- read.delim("mRNA_translation.txt", sep="\t", head=FALSE, col.names = c("gene","description","mRNA_gene")) 
        genes <- read.delim("gene_translation.txt", sep="\t", head=FALSE, col.names = c("gene","description"))

# drop genes with multiple hits depending on isoform (~600 genes)
        combined_list <- left_join(uniprot, mRNAs) %>%
                mutate(gene_and_hit = paste(gene,uniprot_hit, sep = "_"))
        # unique(combined_list$mRNA_gene) %>% length()
        # unique(combined_list$uniprot_hit) %>% length()
        # unique(combined_list$gene) %>% length()
        # unique(combined_list$gene_and_hit) %>% length()
        combined_list2 <- group_by(combined_list, gene) %>%
                summarise(count = n(), hits = length(unique(uniprot_hit))) %>%
                filter(hits == 1)
        unique_hits <- data.frame(gene_hit = unique(combined_list$gene_and_hit)) %>% 
                separate(gene_hit, into = c("gene", "uniprot_hit"))
        final_list <- left_join(combined_list2, unique_hits) %>%
                select(gene, uniprot_hit) %>%
                left_join(genes) %>%
                filter(description != "Gene")
        write.table(final_list, "mya_genes_uniprot_hits.txt", sep="\t", row.names = FALSE, quote = FALSE)

###############
# GENE ONTOLOGY
###############

# Load gene set calls for human (C5: ontology gene sets, C6: oncogenic signature gene sets)
        print("Loading human gene sets")
        C5_gene_set_full <- msigdbr(species = "Homo sapiens", category = "C5") 
        C5_gene_set_full2 <- C5_gene_set_full %>%
                mutate(human_hit = gene_symbol) %>%
                select(human_hit, gs_name) 

# blank tibble
        mya_ncbi_GO <- C5_gene_set_full2[1,] %>%
                mutate(mya_hit = "test") %>%
                filter(mya_hit != "test")

# Loop to make new file (not an efficient way to do this...)
        print("Merging human genes sets with mya arenaria genes - this will take a few minutes")
        total_count <- nrow(final_list)
        count <- 0
        for(g in as.character(final_list$gene)){
                hg <- filter(final_list, gene == g)[,"uniprot_hit"] %>% as.character()
                mya_GO_new <- filter(C5_gene_set_full2, human_hit == hg) %>%
                        mutate(mya_hit = g)
                mya_ncbi_GO <- rbind(mya_ncbi_GO, mya_GO_new)
                count <- count+1
                # print(paste0(count,"/",total_count))
        }
        save(mya_ncbi_GO, file = "mya_ncbi_GO.Rda")

# Make list for running fgsea
        gene_set_summary <- group_by(mya_ncbi_GO[,2:3], gs_name) %>%
                summarise(mya_gene_list = paste0(mya_hit, collapse = " ")) 
        gene_set_list <- as.list(gene_set_summary$mya_gene_list) %>%
                lapply(function(x) unique(strsplit(x, split=" +")[[1]]))
        names(gene_set_list) <- gene_set_summary$gs_name
        save(gene_set_list, file = "mya_gene_set_list.Rda")