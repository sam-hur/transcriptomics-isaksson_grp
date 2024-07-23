#!/bin/bash -l

#SBATCH -A naiss2024-22-952
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 10:00:00
#SBATCH -J Mapping_hisat2
#SBATCH --output=/proj/naiss2024-23-424/analysis/scripts/logs/slurm/slurm-Mapping_hisat2.out
#SBATCH --error=/proj/naiss2024-23-424/analysis/scripts/logs/slurm/slurm-Mapping_hisat2.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user sam.hurenkamp.9631@student.uu.se

source /proj/naiss2024-23-424/analysis/scripts/config.cfg
source /proj/naiss2024-23-424/analysis/scripts/functions.sh

#modules
module load \
    bioinfo-tools \
    HISAT2/2.2.1





