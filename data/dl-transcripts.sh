#!/bin/bash

mkdir -p /proj/data/raw
OUT=/proj/data
SRX_FILE=$OUT/SRX.txt

# Initialize arrays to store accessions
SRX_ACCESSIONS=()
SRR_ACCESSIONS=()

# Read SRX accessions from file
while IFS= read -r line || [ -n "$line" ]; do
    SRX_ACCESSIONS+=("$line")
done < "$SRX_FILE"

# Fetch SRR accessions and download data using prefetch
for SRX_ACCESSION in "${SRX_ACCESSIONS[@]}"; do
    SRR_ACCESSION=$(esearch -db sra -query "$SRX_ACCESSION" | efetch -format runinfo | cut -d ',' -f 1 | grep SRR)
    SRR_ACCESSIONS+=($SRR_ACCESSION)
    for SRR in $SRR_ACCESSION; do
        echo "Prefetching $SRR..."
        prefetch --output-directory $OUT $SRR &
    done
    break
done

# Wait for all prefetch processes to complete
wait

# Log SRR accessions
echo "SRR Accessions: ${SRR_ACCESSIONS[@]}"

# for SRR in "${SRR_ACCESSIONS[@]}"; do
#     if ls "${SRR}"_*.fastq 1> /dev/null 2>&1; then
#         echo "${SRR}_*.fastq already exists, skipping download."
#     else
#         echo "Downloading $SRR..."
#         fasterq-dump --split-files $SRR
#         echo "Zipping $SRR..."
#         gzip ${SRR}_*.fastq
#         rm ${SRR}_*.fastq
#         echo "Deleting $SRR/$SRR.sra"
#         rm -r "$SRR/$SRR.sra"
#         mv "$SRR"*.fastq.gz $SRR
#     fi
# done

find "." -type f -name "*.fastq" | while read -r file; do
    echo "Compressing $file..."
    gzip $file
    # folder=$(basename $(dirname $file))
    # mv *.fastq.gz $folder
done


# Wait for all fasterq-dump processes to complete
wait 

echo "All downloads complete."
