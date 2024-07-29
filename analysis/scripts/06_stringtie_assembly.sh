#!/bin/bash -l

#SBATCH -A naiss2024-22-952
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 16:00:00
#SBATCH -J assembly_stringtie
#SBATCH --output=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-stringtie_assembly.out
#SBATCH --error=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-stringtie_assembly.err
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
    data=$out/hisat2-mapping/$tm
    for sample in $(ls $data); do
        mkdir -p $out/stringtie-assembly/$tm
        
        stringtie \
            -p $T \
            -G $ref_gff \
            -o $out/stringtie-assembly/$tm/merged_transcripts-$sample.$tm.gtf \
            -l $sample \
            $data/$sample/$sample.sorted.bam
    done
done
