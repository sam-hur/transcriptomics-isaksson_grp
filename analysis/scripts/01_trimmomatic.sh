#!/bin/bash -l

#SBATCH -A naiss2024-22-952
#SBATCH -M rackham
#SBATCH -p core
#SBATCH -n 18
#SBATCH -t 12:45:00
#SBATCH -J TRIM
#SBATCH --output=/proj/naiss2024-23-424/analysis/scripts/logs/slurm/slurm-trimmomatic.out
#SBATCH --error=/proj/naiss2024-23-424/analysis/scripts/logs/slurm/slurm-trimmomatic.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user sam.hurenkamp.9631@student.uu.se
source /proj/naiss2024-23-424/transcriptomics-isaksson_grp/analysis/scripts/config.cfg


#modules
module load \
    bioinfo-tools \
    trimmomatic

t=18

for tm in $treatments; do
    echo "processing: $tm"
    
    data_in=$in/raw/$tm
    o=$in/trimmed/$tm
    mkdir -p $o
    
    for dir in $(ls $data_in); do
        echo "-- directory: $dir"
        mkdir -p $o/$dir

        files=($(ls $data_in/$dir))
        raw_R1=$data_in/$dir/${files[0]}
        raw_R2=$data_in/$dir/${files[1]}
        R1P="$o/$dir/${tm}_${dir}_1P.fq.gz"
        R1U="$o/$dir/${tm}_${dir}_1U.fq.gz"
        R2P="$o/$dir/${tm}_${dir}_2P.fq.gz"
        R2U="$o/$dir/${tm}_${dir}_2U.fq.gz"
        
        java -jar $TRIMMOMATIC_ROOT/trimmomatic.jar \
        PE -threads $t -phred33 \
        $raw_R1 $raw_R2 \
        $R1P $R1U $R2P $R2U \
        ILLUMINACLIP:"$TRIMMOMATIC_HOME/adapters/TruSeq3-PE.fa":2:30:7 \
        SLIDINGWINDOW:5:20
    done
done
