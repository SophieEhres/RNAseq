#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --mail-user=${JOB_MAIL}
#SBATCH --mail-type=FAIL
#SBATCH --account=def-sauvagm
#SBATCH --job-name=star_index_firstpass
#SBATCH --ntasks=20
#SBATCH --output=%x_%j.out
#SBATCH --mem-per-cpu=3000

genome=${SCRATCH}/genomes/human/hg38/fasta
index=${SCRATCH}/genomes/human/hg38/STAR_index
fastafile=$(ls ${genome}/*.fa | tr '\n' ' ')

mkdir -p ${index}

STAR --runMode genomeGenerate \
--genomeDir ${index} \
--runThreadN 20 \
--genomeFastaFiles ${fastafile} \

