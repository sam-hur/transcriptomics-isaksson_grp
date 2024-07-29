#!/bin/bash -l

#SBATCH -A naiss2024-22-952
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 12:00:00
#SBATCH -J mapping_hisat2
#SBATCH --output=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-hisat2_mapping.out
#SBATCH --error=/proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/logs/slurm/slurm-hisat2_mapping.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user sam.hurenkamp.9631@student.uu.se

source /proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/config.cfg

#modules
module load \
    bioinfo-tools \
    HISAT2/2.2.1 \
    samtools


# Stop executing pipe on error
set -euo pipefail

T=16

zf_index=$out/hisat2-indexes/Taeniopygia_guttata.index

for tm in $treatments; do
    data=$in/raw/$tm
    for sample in $(ls $data); do
        output=$out/hisat2-mapping/$tm/$sample
        dir=$data/$sample
        mkdir -p $output 
        hisat2 -p $T --dta -x $zf_index \
        -1 $dir/${sample}_1.fastq.gz -2 $dir/${sample}_2.fastq.gz \
        -S $SNIC_TMP/$sample.sam \
        --rna-strandness FR \
        --summary-file $output/$sample-summary.txt

        threads=$((T/2))  # 16 threads causes more overhead, half is more balanced
        samtools view -@ $threads -bS $SNIC_TMP/$sample.sam -bo $SNIC_TMP/$sample.bam
        samtools sort -@ $threads $SNIC_TMP/$sample.bam -o $output/$sample.sorted.bam
        samtools index -@ $threads $output/$sample.sorted.bam -o $output/$sample.sorted.bam.bai
    done
done
