#!/bin/bash -l

#SBATCH -A naiss2024-22-952
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 01:30:00
#SBATCH -J construct_index-hisat2
#SBATCH --output=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-hisat2_build_index.out
#SBATCH --error=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-hisat2_build_index.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user sam.hurenkamp.9631@student.uu.se

source /proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/config.cfg

#modules
module load \
    bioinfo-tools \
    HISAT2/2.2.1

T=16
mkdir -p $out/hisat2-indexes
hisat2-build -p $T $ref $out/hisat2-indexes/Taeniopygia_guttata.index



