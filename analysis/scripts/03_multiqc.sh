#!/bin/bash -l

#SBATCH -A naiss2024-22-952
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -J MULTIQC
#SBATCH --output=/proj/naiss2024-23-424/analysis/scripts/logs/slurm/slurm-MULTIQC.out
#SBATCH --error=/proj/naiss2024-23-424/analysis/scripts/logs/slurm/slurm-MULTIQC.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user sam.hurenkamp.9631@student.uu.se

source /proj/naiss2024-23-424/analysis/scripts/config.cfg
source /proj/naiss2024-23-424/analysis/scripts/functions.sh

#modules
module load \
    bioinfo-tools \
    MultiQC

t=10

i=$out/FastQC  # in
o=$out/MultiQC # out
mkdir -p $o

multiqc $i -o $o -p -v -f

