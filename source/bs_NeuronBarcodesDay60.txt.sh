#!/bin/bash
#SBATCH -n 20
#SBATCH -N 1
#SBATCH -A lsens2018-3-3
#SBATCH -p dell
#SBATCH -t 20:00:00
#SBATCH -J bsplit_NeuronBarcodesDay60.txt
#SBATCH -o bsplit_NeuronBarcodesDay60.txt.out
#SBATCH -e bsplit_NeuronBarcodesDay60.txt.err
# create splitted bam file
module purge
module load subset-bam/1.0

subset-bam -b="/path to your cell barcodes" -o="/path to your bam/possorted_genome_bam.bam"
