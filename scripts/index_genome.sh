#!/bin/bash

# unzip genome and annotation files
    zcat $INPUT_FOLDER/GCF_026914265.1_ASM2691426v1/GCF_026914265.1_ASM2691426v1_genomic.fna.gz \
        > $INPUT_FOLDER/GCF_026914265.1_ASM2691426v1/GCF_026914265.1_ASM2691426v1_genomic.fna
    zcat $INPUT_FOLDER/GCF_026914265.1_ASM2691426v1/GCF_026914265.1_ASM2691426v1_genomic.gtf.gz \
        > $INPUT_FOLDER/GCF_026914265.1_ASM2691426v1/GCF_026914265.1_ASM2691426v1_genomic.gtf

# Index genome for STAR
    # module load star
    STAR \
        --runThreadN 20 \
        --runMode genomeGenerate \
        --genomeDir $INPUT_FOLDER/STARindex \
        --readFilesCommand zcat \
        --genomeFastaFiles $INPUT_FOLDER/GCF_026914265.1_ASM2691426v1/GCF_026914265.1_ASM2691426v1_genomic.fna \
        --sjdbGTFfile $INPUT_FOLDER/GCF_026914265.1_ASM2691426v1/GCF_026914265.1_ASM2691426v1_genomic.gtf