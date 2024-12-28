#!/usr/bin/Rscript

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
#    load("UvsP_BTN_results.rda") # UvsP_BTN: negative is higher in PEI
    load("UvsP_heme_results.rda") # UvsP_heme: negative is higher in PEI
    load("heme_vs_aswheme_results.rda") # heme_ASWheme: negative is higher in untreated

# Load NCBI descriptions and annotations
    gene_location <- paste0(input_folder,"/GCF_026914265.1_ASM2691426v1/gene_translation.txt")
    gene_translations <- read.delim(gene_location,
                                    sep = "\t",
                                    head = FALSE,
                                    col.names = c("gene","description"))
    # uniprot_location <- "/ssd2/Illumina_data/RNAseq/GCF_026914265.1_ASM2691426v1/uniprot_blast/mya_genes_uniprot_hits.txt"
    uniprot_location <- paste0(input_folder,"/GCF_026914265.1_ASM2691426v1/mya_genes_uniprot_hits.txt")
    uniprot_translations <- read.delim(uniprot_location,
                                       sep = "\t",
                                       head = FALSE,
                                       col.names = c("gene", "uniprot", "description")) %>%
        dplyr::select(-description)

# Load gene set annoation information
    # load("/ssd2/Illumina_data/RNAseq/GCF_026914265.1_ASM2691426v1/uniprot_blast/mya_gene_set_list.Rda")
    load(paste0(input_folder,"/GCF_026914265.1_ASM2691426v1/mya_gene_set_list.Rda"))
    C5_gene_set_full <- msigdbr(species = "Homo sapiens", category = "C5") #C5: ontology gene sets
    # load("/ssd2/Illumina_data/RNAseq/GCF_026914265.1_ASM2691426v1/uniprot_blast/mya_ncbi_GO.Rda")
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
        filter(input_annotated, stat < 0, padj < 0.05) %>% #/nrow(input_annotated)) %>%
            arrange(padj) %>%
#            head(n=100) %>%
            write.table(paste0("top100_genes_upreg_in_",comparison,".tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)
        filter(input_annotated, stat > 0, padj < 0.05) %>% #/nrow(input_annotated)) %>%
            arrange(padj)%>%
#            head(n=100) %>%
            write.table(paste0("top100_genes_dnreg_in_",comparison,".tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)
    }

    print_top100(heme_BTN, "BTN_vs_heme")
    print_top100(notheme_BTN, "BTN_vs_tissue")
    print_top100(ASW_BTN, "ASW_vs_BTN")
    print_top100(heme_ASW_BTN, "ASW_vs_heme")
#    print_top100(UvsP_BTN, "USA_vs_PEI_BTN")
    print_top100(UvsP_heme, "USA_vs_PEI_heme")
    print_top100(heme_ASWheme, "heme_vs_ASWheme")

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
                                       stats = ranking)
                                       #nproc = 1) # Changed to try and address issue running for heme_BTN
        gsea_result <- left_join(gsea_result, gs_to_go) %>%
		filter(!str_detect(go_code, "^HP"))
    # print top100 lists
        filter(gsea_result, NES > 0, padj < 0.05) %>%
            arrange(padj) %>%
#            head(n=100) %>%
            dplyr::select(-leadingEdge) %>%
            write.table(paste0("top100_genesets_upreg_in_",comparison,".tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)
        filter(gsea_result, NES < 0, padj < 0.05) %>%
            arrange(padj)%>%
#            head(n=100) %>%
            dplyr::select(-leadingEdge) %>%
            write.table(paste0("top100_genesets_dnreg_in_",comparison,".tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)
        return(gsea_result)
    }
    notheme_BTN_GSEA <- rank_and_GSEA(notheme_BTN, "BTN_vs_tissue")
    save(notheme_BTN_GSEA, file = "notheme_BTN_GSEA.rda")    
    ASW_BTN_GSEA <- rank_and_GSEA(ASW_BTN, "ASW_vs_BTN")
    save(ASW_BTN_GSEA, file = "ASW_BTN_GSEA.rda")
    heme_ASW_BTN_GSEA <- rank_and_GSEA(heme_ASW_BTN, "ASW_vs_heme")
    save(heme_ASW_BTN_GSEA, file = "heme_ASW_BTN_GSEA.rda")
    heme_BTN_GSEA <- rank_and_GSEA(heme_BTN, "BTN_vs_heme")
    save(heme_BTN_GSEA, file = "heme_BTN_GSEA.rda")

#    USA_vs_PEI_BTN_GSEA <- rank_and_GSEA(UvsP_BTN, "USA_vs_PEI_BTN")
#    save(USA_vs_PEI_BTN_GSEA, file = "USA_vs_PEI_BTN_GSEA.rda")
    USA_vs_PEI_heme_GSEA <- rank_and_GSEA(UvsP_heme, "USA_vs_PEI_heme")
    save(USA_vs_PEI_heme_GSEA, file = "USA_vs_PEI_heme_GSEA.rda")
    heme_vs_ASWheme_GSEA <- rank_and_GSEA(heme_ASWheme, "heme_vs_ASWheme")
    save(heme_vs_ASWheme_GSEA, file = "heme_vs_ASWheme_GSEA.rda")

    
    # If GSEA enrichment plots are desired use plotEnrichment function:
    # https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

#### NEW ADDED June 2024 ####
# filter ASW result for only BTN-specific hits    
	asw_heme_notsig <- heme_ASWheme %>%
		as.data.frame() %>%
		filter(!(padj < 0.05)) %>%
		rownames() 
	asw_heme_notsig_BTN <- as.data.frame(ASW_BTN)[asw_heme_notsig,]
    print_top100(asw_heme_notsig_BTN, "ASW_vs_BTN_minus_heme_hits")

# same for GSEA
	asw_heme_notsig_GSEA <- heme_vs_ASWheme_GSEA %>%
		filter(!(padj < 0.05)) %>%
		pull(pathway)
	asw_heme_notsig_BTN_GSEA <- ASW_BTN_GSEA %>%
		filter(pathway %in% asw_heme_notsig_GSEA)

gsea_result <- asw_heme_notsig_BTN_GSEA
comparison <- "ASW_vs_BTN_minus_heme_hits"
        filter(gsea_result, NES > 0, padj < 0.05) %>%
            arrange(padj) %>%
            dplyr::select(-leadingEdge) %>%
            write.table(paste0("top100_genesets_upreg_in_",comparison,".tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)
        filter(gsea_result, NES < 0, padj < 0.05) %>%
            arrange(padj)%>%
            dplyr::select(-leadingEdge) %>%
            write.table(paste0("top100_genesets_dnreg_in_",comparison,".tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)




# Print out top LFC hits
    top_genes <- as.data.frame(heme_BTN) %>%
        rownames_to_column(var = "gene") %>%
        left_join(gene_translations) %>%
        left_join(uniprot_translations) %>%
        dplyr::select(gene, uniprot, padj, baseMean,
                        stat, log2FoldChange, description) %>%
        arrange(padj)

    top_genes %>%
        arrange(log2FoldChange) %>%
            head(n=100) %>%
            write.table(paste0(output_folder,"/plots/top100_LFC_genes_upreg_in_BTN_vs_heme.tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)
    top_genes %>%
        arrange(-log2FoldChange) %>%
            head(n=100) %>%
            write.table(paste0(output_folder,"/plots/top100_LFC_genes_dnreg_in_BTN_vs_heme.tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)

# For supp tables - merge BTNxASW and hemexASW
    # genes
        heme_asw <- as.data.frame(heme_ASWheme) %>%
            rownames_to_column(var = "gene") %>%
            mutate(heme_LFC = -log2FoldChange) %>% # since heme-hemeASW comparison is flipped
            select(gene, heme_LFC, heme_padj = padj)
        asw_treatments <- as.data.frame(ASW_BTN) %>%
            rownames_to_column(var = "gene") %>%
            left_join(gene_translations) %>%
            left_join(uniprot_translations) %>%
            left_join(heme_asw) %>%
            dplyr::select(gene, uniprot, description, baseMean, stat, padj, 
                          log2FoldChange, heme_LFC, heme_padj)
        # print top100 lists
        setwd("/ssd3/RNAseq/outputs/data_tables")
        filter(asw_treatments, stat < 0, padj < 0.05) %>% #/nrow(input_annotated)) %>%
            arrange(padj) %>%
            head(n=100) %>%
            write.table("top100_genes_upreg_in_ASW_vs_BTN_plusheme.tsv",
                        sep="\t", row.names = FALSE, quote = FALSE)
        filter(asw_treatments, stat > 0, padj < 0.05) %>% #/nrow(input_annotated)) %>%
            arrange(padj)%>%
            head(n=100) %>%
            write.table("top100_genes_dnreg_in_ASW_vs_BTN_plusheme.tsv",
                        sep="\t", row.names = FALSE, quote = FALSE)
        asw_treatments %>%
            arrange(padj)%>%
            write.table("ASW_treatment_all_genes_results.tsv",
                        sep="\t", row.names = FALSE, quote = FALSE)

    #gene sets
    heme_asw2 <- heme_vs_ASWheme_GSEA %>%
            mutate(heme_ES = -ES) %>% # since heme-hemeASW comparison is flipped
            dplyr::select(pathway, heme_ES, heme_padj = padj) 

    asw_output <- left_join(ASW_BTN_GSEA, heme_asw2) %>% 
            arrange(padj) %>%
            dplyr::select(-leadingEdge)

        filter(asw_output, NES > 0, padj < 0.05) %>%
            arrange(padj) %>%
            head(n=100) %>%
            write.table("top100_genesets_upreg_in_ASW_vs_BTN_plusheme.tsv",
                        sep="\t", row.names = FALSE, quote = FALSE)
        filter(asw_output,NES < 0, padj < 0.05) %>%
            arrange(padj)%>%
            head(n=100) %>%
            write.table("top100_genesets_dnreg_in_ASW_vs_BTN_plusheme.tsv",
                        sep="\t", row.names = FALSE, quote = FALSE)
        asw_output %>%
            arrange(padj)%>%
            write.table("ASW_treatment_all_genesets_results.tsv",
                        sep="\t", row.names = FALSE, quote = FALSE)


# Compare to mussel results and print top 100
    mussel <- read.delim("./inputs/burioli_mussel_genes.txt",
                         sep = "\t",
                         head = TRUE) %>%
        dplyr::select(uniprot = Gene.ID, Description, log2FC, FDR)

    mussel2 <- inner_join(top_genes, mussel) %>%
        mutate(mya_lfc = -log2FoldChange,
               mya_p_dir = case_when(mya_lfc > 0 ~ -log10(padj), mya_lfc <= 0 ~ log10(padj)),
               mtr_p_dir = case_when(log2FC > 0 ~ -log10(FDR), log2FC <= 0 ~ log10(FDR)))
    mussel3 <- filter(mussel2, padj < 0.05)
    upup <- mussel3 %>% filter(log2FC > 0, mya_lfc > 0) %>% nrow() # 373 both up
    dndn <- mussel3 %>% filter(log2FC < 0, mya_lfc < 0) %>% nrow() # 569 both dn
    updn <- mussel3 %>% filter(log2FC > 0, mya_lfc < 0) %>% nrow() # 286 mussel up, mya dn
    dnup <- mussel3 %>% filter(log2FC < 0, mya_lfc > 0) %>% nrow() # 270 mussel dn, mya up
    mussel3 %>%
        mutate(agg_fold_change = mya_lfc * log2FC) %>%
        arrange(desc(agg_fold_change)) %>%
        head(100) %>%
            write.table(paste0(output_folder,"/plots/top_genes_both_mussels_mya2.tsv"),
                        sep="\t", row.names = FALSE, quote = FALSE)

################# REVISIONS
# compare to inborn errors of immunity

gene_reg_data <- as.data.frame(heme_BTN) %>%
            rownames_to_column(var = "gene") %>%
            left_join(gene_translations) %>%
            left_join(uniprot_translations) %>%
            dplyr::select(gene, uniprot, padj, baseMean,
                          stat, log2FoldChange, description) 
inborn_immunity_genes <- read.delim(paste0(input_folder,"/IUIS-IEI-list-for-web-site-July-2024V2.txt"),
                         sep = "\t",
                         head = TRUE) %>%
		pull(Genetic.defect)
pdf(paste0(output_folder,"/plots/inborn_immunity_genes.pdf"))

dev.off()
