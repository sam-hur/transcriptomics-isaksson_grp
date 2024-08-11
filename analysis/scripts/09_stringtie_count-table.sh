#!/bin/bash -l

#SBATCH -A naiss2024-22-952
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 16:00:00
#SBATCH -J stringtie_counttable
#SBATCH --output=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-stringtie_count_table.out
#SBATCH --error=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-stringtie_count_table.err
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
    ref_gtf=$out/stringtie-merged/TaeGut_transcripts-merged.gtf
    src=$out/hisat2-mapping/$tm
    output=$out/stringtie-count_table/$tm
    mkdir -p $output
    for sample in $(ls $src); do
        echo "processing: $tm, sample: $sample"
        stringtie \
        -p $T \
        -eB \
        -G $ref_gtf \
        -o $output/${sample}.${tm}.count-table.gtf \
        $src/$sample/$sample.sorted.bam
    done
done
