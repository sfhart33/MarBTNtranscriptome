#!/bin/bash

# Script converts genome annotation files to simple tables
    GENOME=$1/GCF_026914265.1_ASM2691426v1
    set -e

# gene name to description file
    # cd $INPUT_FOLDER/GCF_026914265.1_ASM2691426v1
    echo Making gene translation tables with meaningful names
    zcat $GENOME/GCF_026914265.1_ASM2691426v1_genomic.gff.gz | \
        awk -F '\t' ' $3 ~ "^gene" {
            split($9, info, ";" );
            split(info[3], name, "=" );
            split(info[4], description, "=" );
            print name[2], description[2]
            }' OFS="\t" > $GENOME/gene_translation.txt

# mRNA sequence to gene name
    zcat $GENOME/GCF_026914265.1_ASM2691426v1_translated_cds.faa.gz | \
        awk -F '\t' ' $0 ~ "^>" {
            split($0, gene1, "gene=");
            split(gene1[2], gene, "]");
            split($0, protein1, "protein=");
            split(protein1[2], protein, "]");
            split($0, protein_id1, "protein_id=");
            split(protein_id1[2], protein_id, "]");
            print gene[1], protein[1], protein_id[1]
            }' OFS="\t" > $GENOME/mRNA_translation.txt    
    

# Download uniprot database and blast mya sequences for most similar human sequence
    echo Downloading uniprot database
    # Adapted from mya genome paper methods
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz -O $GENOME/uniprot_sprot.fasta.gz
    gunzip $GENOME/uniprot_sprot.fasta.gz

    zcat $GENOME/GCF_026914265.1_ASM2691426v1_protein.faa.gz > \
        $GENOME/GCF_026914265.1_ASM2691426v1_protein.faa

    # module load blast+
    echo Blasting mya sequences for most similar human sequences
    makeblastdb -in $GENOME/uniprot_sprot.fasta -out $GENOME/uniprot_sprot -dbtype prot
    blastp -query $GENOME/GCF_026914265.1_ASM2691426v1_protein.faa -db $GENOME/uniprot_sprot -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out $GENOME/GCF_026914265.1_ASM2691426v1_protein_uniprot_blastp -num_threads $THREADS
    awk -F '\t' ' {
        split($2, columns, "_" );
        split(columns[1], hit, "|" );
        print $1, hit[3]} ' OFS="\t" $GENOME/GCF_026914265.1_ASM2691426v1_protein_uniprot_blastp > $GENOME/GCF_026914265.1_ASM2691426v1_closest_gene_name.txt

# Run r script to merge with human gene ontology databases
    # Rscript  ./scripts/gene_ontology.R