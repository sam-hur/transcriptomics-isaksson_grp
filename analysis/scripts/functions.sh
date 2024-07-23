#!/bin/bash -l

source /proj/naiss2024-23-424/analysis/scripts/config.cfg

get_treatment() {
    # Returns the treatment group for the specified file path, or "" if none is found.
    local filepath="$1"
    local treatments="${2:-$treatments}"

    for treatment in $treatments; 
    do
        if [[ "$filepath" == *"$treatment"* ]]; then
            echo "$treatment"
            return
        fi
    done
    
    echo "No treatment found in filepath"
}