#!/bin/bash
#SBATCH -n 20
#SBATCH -N 1
#SBATCH -A lsens2018-3-3
#SBATCH -p dell
#SBATCH -t 30:00:00
#SBATCH -J countHg38
#SBATCH -o countHg38_%j.out
#SBATCH -e countHg38_%j.err

module purge
module use /projects/fs1/common/modules/
module load GCC/4.9.3-2.25 OpenMPI/1.10.2
module load bcl2fastq2/2.18
module load cellranger/3.0.0

#############################################################
# INPUT PARAMS
#############################################################
TRANSCRIPTOME_DIR="/path to your ReferenceGEnome"
PROJECT="FGF8plus_day16"
SAMPLE="Ale1"
#PROJECTNAME="10XscAle"
#############################################################
# OUTPUT PARAMS
#############################################################
CELLRANGER_ID=$PROJECT"_V3"
FASTQ_PATH="/path to fastq sequencing files"

#############################################################
# Running on multiple cores
#############################################################
NUM_CORES=8
MAX_MEM_PER_CORE=10
MAX_MEM=256
CELLRANGER_OPTIONS="--localcores=$NUM_CORES --mempercore=$MAX_MEM_PER_CORE --localmem=$MAX_MEM"

#############################################################
# MAIN COMMAND
#############################################################
CMD="cellranger count $CELLRANGER_OPTIONS --id=$CELLRANGER_ID --transcriptome=$TRANSCRIPTOME_DIR --fastqs=$FASTQ_PATH --sample=$SAMPLE --project=$PROJECT" 
echo "command: $CMD" 
echo "running on: `uname -a`"
eval $CMD






