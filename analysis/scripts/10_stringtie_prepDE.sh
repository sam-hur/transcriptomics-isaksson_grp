#!/bin/bash -l

#SBATCH -A naiss2024-22-952
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 16:00:00
#SBATCH -J stringtie_prepDE
#SBATCH --output=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-stringtie_count_table.out
#SBATCH --error=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-stringtie_count_table.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user sam.hurenkamp.9631@student.uu.se

source /proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/config.cfg

#modules
module load \
    bioinfo-tools \
    StringTie/2.2.1 \
    python3/3.8.7

# Stop executing pipe on error
set -euo pipefail

T=4 

GTF_SRC=$out/stringtie-prepDE
output=$out/stringtie-prepDE/_csv
mkdir -p $output

prepDE.py -v -i $GTF_SRC -g $output/gene_count_matrix.csv -t $output/transcript_count_matrix.csv
