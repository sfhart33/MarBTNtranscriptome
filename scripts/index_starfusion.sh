#!/bin/bash
    # $1 is $INPUT_FOLDER 
    # $2 is $THREADS

# First download singularity image of STAR-Fusion to $INPUT_FOLDER

# build genome library
    cd $1/GCF_026914265.1_ASM2691426v1
    singularity exec -e -B `pwd` ../star-fusion.v1.11.0.simg \
                /usr/local/src/STAR-Fusion/ctat-genome-lib-builder/prep_genome_lib.pl \
                    --genome_fa ref_genome.fa \
                    --gtf ref_annot.gtf \
                    --dfam_db human \
                    --pfam_db current \
                    --CPU $2
