#!/bin/bash

# Download NCBI-annotated Mya arenaria genome
    wget --recursive --no-host-directories --cut-dirs=6 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/026/914/265/GCF_026914265.1_ASM2691426v1/ -P $INPUT_FOLDER/

# Download RNAseq data from SRA using SRA toolkit
        # Download installation instructions: https://github.com/ncbi/sra-tools/wiki
        # helpful notes on usage https://www.reneshbedre.com/blog/ncbi_sra_toolkit.html
    module load sratools # version 2.10.4

    #test download (for timing)
        # cd $INPUT_FOLDER/fastq
        # fasterq-dump SRR23856993
        # gzip SRR23856993_1.fastq &
        # gzip SRR23856993_2.fastq &
        # wait

    # Full download (takes about an hour each  - will take about 40hr)
	# ASW treatments are archived under same accession number as pre-treatment
        cd $INPUT_FOLDER/fastq
        ACCESSIONS="SRR22741477
            SRR22741478
            SRR22741479
            SRR22741480
            SRR22741481
            SRR22741482
            SRR23856977
            SRR23856979
            SRR23856980
            SRR23856981
            SRR23856982
            SRR23856983
            SRR23856984
            SRR23856985
            SRR23856986
            SRR23856987
            SRR23856988
            SRR23856989
            SRR23856990
            SRR23856991
            SRR23856992
            SRR23856993
            SRR23856994
            SRR23856995
            SRR30420898
            SRR30420899
            SRR30420900"
        TOTAL=$(printf '%s\n' $ACCESSIONS:q | wc -w)
        COUNT=0
        for ACCESSION in $ACCESSIONS
        do
            echo Downloading $ACCESSION
            COUNT=$((COUNT+1))
            fasterq-dump $ACCESSION
            echo Zipping $ACCESSION
            gzip $ACCESSION"_1.fastq" &
            gzip $ACCESSION"_2.fastq" &
            echo $COUNT complete out of $TOTAL
        done
        wait