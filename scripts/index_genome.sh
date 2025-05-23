#!/bin/bash
    # $1 is $INPUT_FOLDER 
    # $2 is $THREADS

# unzip genome and annotation files
    # zcat $1/GCF_026914265.1_ASM2691426v1/GCF_026914265.1_ASM2691426v1_genomic.fna.gz \
    #     > $1/GCF_026914265.1_ASM2691426v1/GCF_026914265.1_ASM2691426v1_genomic.fna
    # zcat $1/GCF_026914265.1_ASM2691426v1/GCF_026914265.1_ASM2691426v1_genomic.gtf.gz \
    #     > $1/GCF_026914265.1_ASM2691426v1/GCF_026914265.1_ASM2691426v1_genomic.gtf

    zcat $1/GCF_026914265.1_ASM2691426v1/GCF_026914265.1_ASM2691426v1_genomic.fna.gz \
        > $1/GCF_026914265.1_ASM2691426v1/ref_genome.fa
    zcat $1/GCF_026914265.1_ASM2691426v1/GCF_026914265.1_ASM2691426v1_genomic.gtf.gz \
        > $1/GCF_026914265.1_ASM2691426v1/ref_annot.gtf

# Index genome for STAR
    # module load star
    echo Indexing mya arenaria genome for STAR
    STAR \
        --runThreadN $2 \
        --runMode genomeGenerate \
        --genomeDir $1/STARindex \
        --readFilesCommand zcat \
        --genomeFastaFiles $1/GCF_026914265.1_ASM2691426v1/ref_genome.fa \
        --sjdbGTFfile $1/GCF_026914265.1_ASM2691426v1/ref_annot.gtf