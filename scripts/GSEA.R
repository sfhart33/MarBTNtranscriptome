# Script to run Gene set enrichment analysis

# input variables
    input_folder <- commandArgs(trailingOnly=TRUE)[1]
    output_folder <- commandArgs(trailingOnly=TRUE)[2]

# Load packages
    library(tidyverse)
    library(DESeq2)
    library(fgsea)
    library(msigdbr)
    library(stats)
    library(GO.db)

# Load deseq results objects
    setwd(paste0(output_folder,"/data_tables"))
    load("asw_vs_btn_results.rda") # ASW_BTN: negative is higher in ASW
    load("heme_vs_ASWbtn_results.rda") # heme_BTN: negative is higher in BTN
    load("heme_vs_btn_results.rda") # notheme_BTN: negative is higher in BTN
    load("notheme_vs_btn_results.rda") # heme_ASW_BTN: negative is higher in BTN-ASW

# Load NCBI descriptions and annotations
    gene_location <- paste0(input_folder,"/GCF_026914265.1_ASM2691426v1/gene_translation.txt")
    gene_translations <- read.delim(gene_location,
                                    sep = "\t",
                                    head = FALSE,
                                    col.names = c("gene","description"))
    #uniprot_location <- "/ssd2/Illumina_data/RNAseq/GCF_026914265.1_ASM2691426v1/uniprot_blast/mya_genes_uniprot_hits.txt"
    uniprot_location <- paste0(input_folder,"/GCF_026914265.1_ASM2691426v1/mya_genes_uniprot_hits.txt")
    uniprot_translations <- read.delim(uniprot_location,
                                       sep = "\t",
                                       head = FALSE,
                                       col.names = c("gene", "uniprot", "description")) %>%
        dplyr::select(-description)

# Load gene set annoation information
    load(paste0(input_folder,"/GCF_026914265.1_ASM2691426v1/mya_gene_set_list.Rda"))
    C5_gene_set_full <- msigdbr(species = "Homo sapiens", category = "C5") #C5: ontology gene sets
    load(paste0(input_folder,"/GCF_026914265.1_ASM2691426v1/mya_ncbi_GO.Rda"))
    count_GOs <- mya_ncbi_GO$gs_name %>% unique() %>% length()
    gs_info_df <- C5_gene_set_full %>%
        dplyr::select(gs_subcat,gs_name,
                      gs_exact_source,
                      gs_url,
                      gs_description) %>%
        distinct()
    count_genesets <- nrow(gs_info_df)
    gs_to_go <- data.frame(pathway = unique(C5_gene_set_full$gs_name),
                           go_code = unique(C5_gene_set_full$gs_exact_source))

# Merge annotations and print top 100 lists
    print_top100 <- function(input, comparison){
        input_annotated <- as.data.frame(input) %>%
            rownames_to_column(var = "gene") %>%
            left_join(gene_translations) %>%
            left_join(uniprot_translations) %>%
            dplyr::select(gene, uniprot, padj, baseMean,
                          stat, log2FoldChange, description)
        # print top100 lists
        filter(input_annotated, stat < 0, padj < 0.05/nrow(input_annotated)) %>%
            arrange(padj) %>%
            head(n=100) %>%
            write.table(paste0("top100_genes_upreg_in_",comparison,".tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)
        filter(input_annotated, stat > 0, padj < 0.05/nrow(input_annotated)) %>%
            arrange(padj)%>%
            head(n=100) %>%
            write.table(paste0("top100_genes_dnreg_in_",comparison,".tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)
    }

    print_top100(heme_BTN, "BTN_vs_heme")
    print_top100(notheme_BTN, "BTN_vs_tissue")
    print_top100(ASW_BTN, "ASW_vs_BTN")
    print_top100(heme_ASW_BTN, "ASW_vs_heme")


# Run Gene Set Enrichment Analysis
    rank_and_GSEA <- function(input, comparison){
    # ranks genes from most up to down regulated
        set.seed(12345)
        ranking <- input %>%
            as.data.frame() %>%
            rownames_to_column(var = "gene") %>%
            distinct(gene, .keep_all = TRUE) %>%
            filter(!is.na(stat)) %>%
            mutate(rank = rank(-stat, ties.method = "random")) %>%
            arrange(desc(rank)) %>% 
            mutate(statrank = -stat) %>%
            dplyr::select(gene, statrank) %>%
            #select(gene, rank) %>%
            deframe()
    # Run GSEA
        set.seed(12345)
        gsea_result <- fgseaMultilevel(pathways = gene_set_list,
                        stats    = ranking) %>%
            left_join(gs_to_go)
    # print top100 lists
        filter(gsea_result, NES > 0) %>%
            arrange(padj) %>%
            head(n=100) %>%
            dplyr::select(-leadingEdge) %>%
            write.table(paste0("top100_genesets_upreg_in_",comparison,".tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)
        filter(gsea_result, NES < 0) %>%
            arrange(padj)%>%
            head(n=100) %>%
            dplyr::select(-leadingEdge) %>%
            write.table(paste0("top100_genesets_dnreg_in_",comparison,".tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)
        return(gsea_result)
    }
    ASW_BTN_GSEA <- rank_and_GSEA(ASW_BTN, "ASW_vs_BTN")
    save(ASW_BTN_GSEA, file = "ASW_BTN_GSEA.rda")
    heme_ASW_BTN_GSEA <- rank_and_GSEA(heme_ASW_BTN, "ASW_vs_heme")
    save(heme_ASW_BTN_GSEA, file = "heme_ASW_BTN_GSEA.rda")
    heme_BTN_GSEA <- rank_and_GSEA(heme_BTN, "BTN_vs_heme")
    save(heme_BTN_GSEA, file = "heme_BTN_GSEA.rda")
    notheme_BTN_GSEA <- rank_and_GSEA(notheme_BTN, "BTN_vs_tissue")
    save(notheme_BTN_GSEA, file = "notheme_BTN_GSEA.rda")

    # If GSEA enrichment plots are desired use plotEnrichment function:
    # https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html


