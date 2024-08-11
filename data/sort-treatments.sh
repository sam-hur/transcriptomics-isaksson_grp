#!/bin/bash

content=$(cat runs.txt)
treatments="ALAN soot noise control"

mkdir -p $treatments

echo "$content" | while IFS=',' read -r run_id treatment srr_code
do
    if [ -d "$srr_code" ]; then
        echo "Processing: $srr_code --> $treatment"
        mv "$srr_code" "$treatment/"
        echo "Moved $srr_code to $treatment/"
    else
        echo "Directory $srr_code does not exist."
    fi
done


for tm in $treatments; do
    data_in=$in/raw/$tm
    find $data_in -mindepth 2 -type f -exec mv '{}' $data_in ';'
    find $data_in -type d -empty -delete
done
