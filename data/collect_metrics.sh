#!/bin/bash

# Initialize an empty array
echo "[]" > complete_metrics.json

# Loop through each directory
for d in */*/metrics; do
    path=$d

    # Extract 'tissue' and 'obs' values from the path
    tissue=$(echo $path | cut -d'/' -f1)
    obs=$(echo $path | cut -d'/' -f2)

    # Combine JSON files
    jq -s '.[0] * .[1] * .[2]' "$d/run_info.json" "$d/raw_metrics.json" "$d/kb_info.json" | jq --arg tissue "$tissue" --arg obs "$obs" '. + {tissue: $tissue, obs: $obs}' > $d/combined.json

done;



# Find all 'combined.json' files and process them
find . -type f -name 'combined.json' | while read -r file; do
    if jq -e 'type == "array"' "$file" > /dev/null; then
        # File is an array, merge its elements
        jq -s '.[0] + .[1]' ./complete_metrics.json "$file" > ./all_tmp.json && mv ./all_tmp.json ./complete_metrics.json
    elif jq -e 'type == "object"' "$file" > /dev/null; then
        # File is an object, add it as an element
        jq '. += [$new]' --argjson new "$(cat "$file")" ./complete_metrics.json > ./all_tmp.json && mv ./all_tmp.json ./complete_metrics.json
    else
        echo "Error: JSON file $file is neither an array nor an object."
    fi
done

mv ./complete_metrics.json metrics.json
