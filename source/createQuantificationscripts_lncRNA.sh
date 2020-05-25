#!/bin/bash
#SBATCH -n 20
#SBATCH -N 1
#SBATCH -A lsens2018-3-3
#SBATCH -p dell
#SBATCH -t 19:00:00
#SBATCH -J FCcount
#SBATCH -o %j.FCcount.out
#SBATCH -e %j.FCcount.err


ml GCC/5.4.0-2.26 OpenMPI/1.10.3 Subread/1.5.0-p1

files=$(ls “/path to your bam/possorted_genome_bam.bam)

# Gencode exon
featureCounts -T 8 -p -B \
  --primary -F GTF -t exon -g gene_name \
    -a “/path to your bam/possorted_Reference genomegenome_bam.bam \
    -o /projects/fs1/yogitas/Seq010_Oct30_RE/180604_NB502004_0018_AHYMM3BGX5/Cellbarcodes/FeatureCounts/SCells_LncRNA_s1.txt $files
