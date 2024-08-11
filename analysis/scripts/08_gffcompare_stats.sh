#!/bin/bash -l

#SBATCH -A naiss2024-22-952
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -J stats-gffcompare
#SBATCH --output=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-gffcompare.out
#SBATCH --error=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-gffcompare.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user sam.hurenkamp.9631@student.uu.se

source /proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/config.cfg

#modules
# gffcompare is installed locally and sits in my bashrc env variables.
gffcompare=/home/samhur/gffcompare/gffcompare

# Stop executing pipe on error
set -euo pipefail

# for tm in $treatments; do
#     data=$out/gffcompare-stats/$tm
#     mkdir -p $data
#     $gffcompare -V -G -r $ref_gff -o $data/merged-$tm.gtf $out/stringtie-merged/$tm/*.gtf
# done
output=$out/gffcompare-stats

# grab total annotations
# files=$(find $out/gffcompare-stats -type f -name "*.gtf")
# $gffcompare -V -G -r $ref_gff -o $output/merged-treatments.gtf $out/stringtie-merged/**/*.gtf

mkdir -p $output/TaeGut-merged
$gffcompare -V -G -r $ref_gff -o $output/TaeGut-merged/merged-TaeGut.gtf $out/stringtie-merged/TaeGut_transcripts-merged.gtf
