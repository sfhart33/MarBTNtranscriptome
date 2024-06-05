#!/usr/bin/Rscript

# Script to analysis fusion gene outputs

# Load packages and set wd
    library(tidyverse)
    library(DESeq2)
    library(gridExtra)

# input variables
    input_folder <- commandArgs(trailingOnly=TRUE)[1]
    output_folder <- commandArgs(trailingOnly=TRUE)[2]
    samples_data <- read.table("./inputs/sample_data.txt", header = TRUE) 
    home_dir <- getwd()
    setwd(paste0(output_folder,"/star_fusion"))

# for testing
    # input_folder <- "/ssd3/RNAseq/inputs"
    # output_folder <- "/ssd3/RNAseq/outputs"
    # samples_data <- read.table("./inputs/sample_data.txt", header = TRUE) 
    # home_dir <- getwd()
    # setwd(paste0(output_folder,"/star_fusion"))

# Set groupings
    samples_names <- filter(samples_data, tissue %in% c("heme","BTN")) %>%
        pull(name) %>%
        as.character()
    samples_clams <- filter(samples_data, tissue %in% c("heme","BTN")) %>%
        pull(clam) %>%
        as.character()
    samples_type <- filter(samples_data, tissue %in% c("heme","BTN")) %>%
        pull(tissue) %>%
        as.character()
    samples_location <- filter(samples_data, tissue %in% c("heme","BTN")) %>%
        pull(location) %>%
        as.character()
    samples_type[samples_type=="heme"] <- "healthy"
    BTN_samples <- filter(samples_data, tissue == "BTN") %>%
        pull(name) %>%
        as.character()
    BTN_clams <- filter(samples_data, tissue == "BTN") %>%
        pull(clam) %>%
        as.character()
    heme_samples <- filter(samples_data, tissue == "heme") %>%
        pull(name) %>%
        as.character()
    heme_clams <- filter(samples_data, tissue == "heme") %>%
        pull(clam) %>%
        as.character()
    USA_BTN_clams <- filter(samples_data, tissue == "BTN", location == "USA") %>%
        pull(clam) %>%
        as.character()    
#    PEI_BTN_clams <- filter(samples_data, tissue == "BTN", location == "PEI")  %>%
#        pull(clam) %>%
#        as.character()

# Load data and count number shared by all cancers
  for (sample in samples_names){
        sample_load <- read.table(
                        # paste0("/ssd3/RNAseq/fusions_ncbi/fusion_analysis/",sample,".fusion_predictions.abridged.tsv"),
                        paste0(sample,".fusion_predictions.abridged.tsv"),
                        header = TRUE,
                        comment.char = "") %>%
            select( LeftGene, RightGene, LeftBreakpoint, RightBreakpoint, FFPM)
        colnames(sample_load) <- c("LeftGene", "RightGene", "LeftBreakpoint", "RightBreakpoint", sample)
        # print(sample)
        # print(nrow(sample_load))
        if (sample == samples_names[1]){
            fusions_file <- sample_load 
        } else {
            fusions_file <- full_join(fusions_file, sample_load)
        }
        # print(nrow(fusions_file))
    }
    # head(fusions_file)
    # nrow(fusions_file)
    fusions_file[is.na(fusions_file)] <- 0
    colnames(fusions_file) <- c("LeftGene", "RightGene", "LeftBreakpoint", "RightBreakpoint", samples_clams)

# Count overlp among samples
    cancer_shared <- fusions_file %>%
        select(all_of(BTN_clams)) %>%
        filter_all(all_vars(. > 0)) %>%
        nrow() # 212 in all BTN

#    usabtn_shared <- fusions_file %>%
#        select(all_of(USA_BTN_clams)) %>%
#        filter_all(all_vars(. > 0)) %>%
#        nrow() # 212 in all USA BTN

#    usapei_shared <- fusions_file %>%
#        filter(BTN_6> 0) %>% # update with other PEI
#        select(all_of(USA_BTN_clams)) %>%
#        filter_all(any_vars(. > 0)) %>%
#        nrow() # 91 in any USA BTN and PEI BTN

    cancer_shared_nohealthy <- fusions_file %>%
        filter_at(vars(starts_with("BTN")), all_vars(. > 0)) %>%
        filter_at(vars(starts_with("HEL")), all_vars(. == 0)) %>%
        nrow() # 181 in all BTN no healthy 

    usabtn_shared_nohealthy <- fusions_file %>%
        filter_at(vars(starts_with("HEL")), all_vars(. == 0)) %>%
        select(all_of(USA_BTN_clams)) %>%
        filter_all(all_vars(. > 0)) %>%
        nrow() # 181 in all USA BTN no healthy

    all_sample_fusions <- fusions_file %>%
        select(all_of(samples_clams)) %>%
        filter_all(all_vars(. > 0)) %>%
        nrow() # 16 in all samples

    all_healthy_fusions <- fusions_file %>%
        select(all_of(heme_clams)) %>%
        filter_all(all_vars(. > 0)) %>%
        nrow() # 16 in all healthies

# Tally up total fusions per sample
    fusion_counts <- fusions_file %>%
        select(samples_clams) %>%
        t() %>%
        as.data.frame() %>%
        mutate(count=rowSums(.!=0)) %>%
        select(count) %>%
        add_column(type = samples_type, location = samples_location) %>%
        mutate(group = ifelse(type == "BTN", paste(location, type), type))

    # fusion_counts_mean <- fusion_counts %>%
    #     group_by(group) %>%
    #     summarize(mean = mean(count),
    #             sd = sd(count))

# Correct for fusions in all samples
    fusion_counts_adj <- mutate(fusion_counts, adj.count = as.numeric(count) - all_sample_fusions)
    fusion_counts_mean_adj <- fusion_counts_adj %>%
        group_by(type) %>%
        summarize(mean = mean(adj.count),
                sd = sd(adj.count))
    #fusion_counts_mean$shared <- c(cancer_shared,all_healthy_fusions)
    fusion_counts_mean_adj$shared <- c(cancer_shared - all_sample_fusions, 0)

# t-test
    t.test(as.numeric(filter(fusion_counts, type=="healthy")[,"count"]),as.numeric(filter(fusion_counts, type=="BTN")[,"count"])) # p-value = 0.0001548
    # t.test(as.numeric(filter(fusion_counts, group=="USA BTN")[,"count"]),as.numeric(filter(fusion_counts, group=="PEI BTN")[,"count"])) # not enough 'y' observations
    

# merge genes with human-readable names
    cancer_shared_nohealthy_df <-  fusions_file %>%
        filter_at(vars(starts_with("BTN")), all_vars(. > 0)) %>%
        filter_at(vars(starts_with("HEL")), all_vars(. == 0))
    gene_names1 <- read.table(paste0(input_folder,"/GCF_026914265.1_ASM2691426v1/gene_translation.txt"),
        sep="\t", header = FALSE, col.names = c("LeftGene","LeftGene_desc"), quote = "")
    gene_names2 <- read.table(paste0(input_folder,"/GCF_026914265.1_ASM2691426v1/gene_translation.txt"),
        sep="\t", header = FALSE, col.names = c("RightGene","RightGene_desc"), quote = "")
    fusion_gene_list <- select(cancer_shared_nohealthy_df, LeftGene, RightGene) %>%
        left_join(gene_names1) %>%
        left_join(gene_names2)
    #fusion_gene_list %>% head()
    write.table(fusion_gene_list, "fusion_genes_list.txt", row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")

# Check for kinases
    filter(fusion_gene_list, grepl("kinase",RightGene_desc) | grepl("kinase",LeftGene_desc)) # 8 kinases out of 258 unique genes involved in fusions = 3.1%
    filter(gene_names2, grepl("kinase",RightGene_desc)) %>% nrow() # 621/43039 = 1.44%
    gene_names2 %>% nrow() 
    # chi quared = 0.025542 (excel)
    # filter(fusion_gene_list, grepl("kinase",RightGene_desc) | grepl("kinase",LeftGene_desc)) %>%
    #     write.table("fusion_kinases_list.txt", row.names = TRUE, col.names = FALSE, quote = FALSE)

# plot
  pdf(paste0(output_folder,"/plots/fusion_plot.pdf"))
        ggplot(fusion_counts_mean_adj, aes(x=type,y=mean)) +
            geom_bar(position="dodge", stat="identity", fill = "#969696")+
            geom_bar(aes(x=type,y=shared),position="dodge", stat="identity", fill = "black")+
            geom_errorbar(aes(ymin=(mean-sd), ymax=(mean+sd)), width=.2,position=position_dodge(.9))+
            geom_point(data=fusion_counts_adj,
                        aes(x=type,y=as.numeric(adj.count)),
                        position=position_dodge(.9), size=2)+
            xlab(NULL)+
            ylab("Number of fusion transcripts\n(minus fusions in all samples)")+
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold"),
                    legend.title = element_blank())
    dev.off()

#########################################################################################################
# CN

setwd(output_folder)

# load gene by cn list
    usa_calls <- paste0(home_dir,"/inputs/genes_by_USA_copy_number.tsv") %>%
        read.delim(header = FALSE, col.names = c("scaf", "start", "end", "mya_gene", "cn")) %>%
        separate(mya_gene, into = c(NA, "gene"), sep = "gene-", ) %>%
        select(gene, cn) %>%
        filter(cn != 0)
    pei_calls <- paste0(home_dir,"/inputs/genes_by_PEI_copy_number.tsv") %>%
        read.delim(header = FALSE, col.names = c("scaf", "start", "end", "mya_gene", "cn")) %>%
        separate(mya_gene, into = c(NA, "gene"), sep = "gene-", ) %>%
        select(gene, cn) %>%
        filter(cn != 0)

# Load gene annotation list
    gene_names <- read.table(paste0(input_folder,"/GCF_026914265.1_ASM2691426v1/gene_translation.txt"),
        sep="\t", header = FALSE, col.names = c("gene","gene_desc"), quote = "") 
    head(gene_names)

# Load expression data and merge
#    load(paste0(output_folder,"/data_tables/heme_vs_usabtn_results.rda")) # heme_uBTN: negative is higher in BTN
#    usa_exp <- as.data.frame(heme_uBTN) %>%
#        rownames_to_column(var = "gene") %>%
#        filter(baseMean > 1) %>%
#        select(gene, log2FoldChange) %>%
#        inner_join(usa_calls) %>%
#        mutate(sublin = "USA")
#    load(paste0(output_folder,"/data_tables/heme_vs_peibtn_results.rda")) # heme_pBTN: negative is higher in BTN    
#    pei_exp <- as.data.frame(heme_pBTN) %>%
#        rownames_to_column(var = "gene") %>%
#        filter(baseMean > 1) %>%
#        select(gene, log2FoldChange) %>%
#        inner_join(pei_calls) %>%
#        mutate(sublin = "PEI")

#    cn_exp <- rbind(usa_exp, pei_exp)

    load(paste0(output_folder,"/data_tables/heme_vs_btn_results.rda")) # heme_BTN: negative is higher in BTN
    cn_exp <- as.data.frame(heme_BTN) %>%
        rownames_to_column(var = "gene") %>%
        filter(baseMean > 1) %>%
        select(gene, log2FoldChange) %>%
        inner_join(usa_calls) %>%
        mutate(sublin = "USA")



pdf(paste0(output_folder,"/plots/cn_exp_plot.pdf"))
    ggplot(cn_exp, aes(cn, -log2FoldChange)) +
        geom_hline(yintercept = 0) +  
        geom_boxplot(outlier.shape = NA) +
        ylim(-12,12) +
        ylab("log2(foldchange) in BTN vs hemocytes")+
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold"),
                    legend.title = element_blank())
dev.off()

#########################################################################################################
# Steamer
    # load("/ssd3/RNAseq/outputs/data_tables/heme_vs_btn_results.rda") # heme_BTN: negative is higher in BTN

# Load steamer sites from genome paper data, join with expression data
    steamer_ins <- paste0(home_dir,"/inputs/steamer_genes.tsv") %>%
        read.delim(header = FALSE, col.names = c("gene", "inserted", "strand", "group")) %>%
        left_join(gene_names) %>%
        left_join(
            dplyr::select(
                rownames_to_column(
                    as.data.frame(heme_BTN), var = "gene"
                ), gene, u_lfc = log2FoldChange, u_padj = padj)
            ) %>%         
        # mutate(
        #     U_up = ifelse(u_lfc < 0, TRUE, FALSE),
        #     U_dn = ifelse(u_lfc > 0, TRUE, FALSE),
        #     P_up = ifelse(p_lfc < 0, TRUE, FALSE),
        #     P_dn = ifelse(p_lfc > 0, TRUE, FALSE))        
        mutate(
            U_up = ifelse(u_padj < 0.05 & u_lfc < 0, TRUE, FALSE),
            U_dn = ifelse(u_padj < 0.05 & u_lfc > 0, TRUE, FALSE),
#            P_up = ifelse(p_padj < 0.05 & p_lfc < 0, TRUE, FALSE),
#            P_dn = ifelse(p_padj < 0.05 & p_lfc > 0, TRUE, FALSE),
            PEI = ifelse(group %in% c("allCnoH","allPEI"), TRUE, FALSE),
            USA = ifelse(group %in% c("allCnoH","allUSA"), TRUE, FALSE))

# Count up/down regulated genes
    steamer_upstream_USA <- steamer_ins %>%
        filter(inserted == "2kbup") %>%
        group_by(strand, USA) %>%
        summarise(
            up=sum(U_up, na.rm = TRUE),
            dn=sum(U_dn, na.rm = TRUE),
#            up_con=sum(P_up, na.rm = TRUE),
#            dn_con=sum(P_dn, na.rm = TRUE),
            n=n()
            )

#    steamer_upstream_PEI <- steamer_ins %>%
#        filter(inserted == "2kbup") %>%
#        group_by(strand, PEI) %>%
#        summarise(
#            up=sum(P_up, na.rm = TRUE),
#            dn=sum(P_dn, na.rm = TRUE),
#            up_con=sum(U_up, na.rm = TRUE),
#            dn_con=sum(U_dn, na.rm = TRUE),
#            n=n()
#            )

    steamer_genes_USA <- steamer_ins %>%
        filter(inserted == "genes") %>%
        group_by(USA) %>%
        summarise(
            up=sum(U_up, na.rm = TRUE),
            dn=sum(U_dn, na.rm = TRUE),
#            up_con=sum(P_up, na.rm = TRUE),
#            dn_con=sum(P_dn, na.rm = TRUE),
            n=n()
            )

#    steamer_genes_PEI <- steamer_ins %>%
#        filter(inserted == "genes") %>%
#        group_by(PEI) %>%
#        summarise(
#            up=sum(P_up, na.rm = TRUE),
#            dn=sum(P_dn, na.rm = TRUE),
#            up_con=sum(U_up, na.rm = TRUE),
#            dn_con=sum(U_dn, na.rm = TRUE),
#            n=n()
#            )

# steamer_upstream_USA
# steamer_upstream_PEI
# steamer_genes_USA
# steamer_genes_PEI

#### NEW boxplot of LFC based on steamer insertion status

pdf(paste0(output_folder,"/plots/steamer_exp_plot.pdf"))
    ggplot(steamer_ins, aes(inserted, -u_lfc, fill = USA)) +
        geom_hline(yintercept = 0) +  
        geom_boxplot(outlier.shape = NA) +
        #ylim(-12,12) +
        ylab("log2(foldchange) in BTN vs hemocytes")+
	xlab("Present in samples")+
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold"),
                    legend.title = element_blank())
    ggplot(steamer_ins, aes(inserted, -u_lfc, fill = USA, color = strand)) +
        geom_hline(yintercept = 0) +  
        geom_boxplot(outlier.shape = NA) +
	scale_fill_manual(values = c("green","orange")) +
        #ylim(-12,12) +
        ylab("log2(foldchange) in BTN vs hemocytes")+
	xlab("Present in samples")+
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold"),
                    legend.title = element_blank())

    ggplot(steamer_ins, aes(inserted, -u_lfc, fill = USA)) +
        geom_hline(yintercept = 0) +  
        geom_violin() +
        #ylim(-12,12)
        ylab("log2(foldchange) in BTN vs hemocytes")+
	xlab("Present in samples")+
            theme_classic() +
            theme(axis.text=element_text(size=12,face="bold"),
                    axis.title=element_text(size=16,face="bold"),
                    text=element_text(size=18,face="bold"),
                    legend.title = element_blank())

dev.off()

t.test(filter(steamer_ins, inserted=="2kbup", USA==TRUE)[,"u_lfc"],filter(steamer_ins, inserted=="2kbup", USA==FALSE)[,"u_lfc"]) # p-value = 0.1623
t.test(filter(steamer_ins, inserted=="2kbup", USA==TRUE)[,"u_lfc"]) # p-value = 0.01015
t.test(filter(steamer_ins, inserted=="2kbup", USA==FALSE)[,"u_lfc"]) # p-value = 0.5117
wilcox.test(filter(steamer_ins, inserted=="2kbup", USA==TRUE)[,"u_lfc"],filter(steamer_ins, inserted=="2kbup", USA==FALSE)[,"u_lfc"]) # p-value = 0.0384

t.test(filter(steamer_ins, inserted=="genes", USA==TRUE)[,"u_lfc"],filter(steamer_ins, inserted=="genes", USA==FALSE)[,"u_lfc"]) # p-value = 0.385
t.test(filter(steamer_ins, inserted=="genes", USA==TRUE)[,"u_lfc"]) # p-value = 0.02075
t.test(filter(steamer_ins, inserted=="genes", USA==FALSE)[,"u_lfc"]) # p-value = 0.06387
wilcox.test(filter(steamer_ins, inserted=="genes", USA==TRUE)[,"u_lfc"],filter(steamer_ins, inserted=="genes", USA==FALSE)[,"u_lfc"]) # p-value = 0.2679

# Print output tables
    # Did not plot since there doesn't appear to be a significant effect of steamer with the new samples
#    cbind(steamer_genes_USA[,1],
#        steamer_genes_USA[,-1] + steamer_genes_PEI[,-1]) %>%
#        mutate(
#            up_frac = up/(up+dn),
#            up_frac_con = up_con/(up_con+dn_con),
#            up_frac2 = up/n,
#            up_frac_con2 = up_con/n
#        ) %>%
#        write.table(paste0(output_folder,"/data_tables/steamer_ins_in_gene_effect.tsv"),
#                    sep="\t", row.names = TRUE, quote = FALSE)
#    cbind(ungroup(steamer_upstream_USA)[,1:2],
#        steamer_upstream_USA[,-c(1,2)] + steamer_upstream_PEI[,-c(1,2)]) %>%
#        mutate(
#            up_frac = up/(up+dn),
#            up_frac_con = up_con/(up_con+dn_con),
#            up_frac2 = up/n,
#            up_frac_con2 = up_con/n
#        ) %>%
#        write.table(paste0(output_folder,"/data_tables/steamer_ins_upstream_gene_effect.tsv"),
#                    sep="\t", row.names = TRUE, quote = FALSE)