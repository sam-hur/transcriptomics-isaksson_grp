#!/bin/bash -l

#SBATCH -A naiss2024-22-952
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 01:45:00
#SBATCH -J QC_POST-TRIM
#SBATCH --output=/proj/naiss2024-23-424/analysis/scripts/logs/slurm/slurm-fastqc_post.out
#SBATCH --error=/proj/naiss2024-23-424/analysis/scripts/logs/slurm/slurm-fastqc_post.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user sam.hurenkamp.9631@student.uu.se
source /proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/config.cfg

#modules
module load \
    bioinfo-tools \
    FastQC

t=10

for tm in $treatments; do
    echo "processing: $tm"
    o=$out/FastQC/post_trim/$tm
    mkdir -p $o
    fastqc -o $o -t $t $in/trimmed/$tm/*/*.fq.gz
done
