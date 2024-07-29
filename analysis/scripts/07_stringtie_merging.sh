#!/bin/bash -l

#SBATCH -A naiss2024-22-952
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 16:00:00
#SBATCH -J merging_stringtie
#SBATCH --output=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-stringtie_merging.out
#SBATCH --error=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-stringtie_merging.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user sam.hurenkamp.9631@student.uu.se

source /proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/config.cfg

#modules
module load \
    bioinfo-tools \
    StringTie/2.2.1 \


# Stop executing pipe on error
set -euo pipefail

T=16

for tm in $treatments; do
    data=$out/stringtie-assembly/$tm
    
    mkdir -p $out/stringtie-merged/$tm
    
    stringtie \
        --merge -p $T \
        -G $ref_gff \
        -o $out/stringtie-merged/$tm/$tm-merged.gtf \
        $(ls $data/*.gtf)
done
