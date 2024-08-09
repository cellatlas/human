#!/bin/bash

# Function to extract links from observation.json
extract_links() {
    local json_file="$1"
    jq -r '.links[] | "\(.url) \(.filetype) \(.accession) \(.filename)"' "$json_file"
}

# Function to process each link
process_link() {
    local url="$1"
    local filetype="$2"
    local accession="$3"
    local filename="$4"
    local output_dir="$5"

    mkdir -p "$output_dir"

    case "$filetype" in
        "bam")
            local filename=$(basename "$url")
            wget --continue --quiet -O "$output_dir/$filename" "$url"
            bamtofastq --reads-per-fastq=500000000 "$output_dir/$filename" "$output_dir/fastq"
            rm "$output_dir/$filename"  # Remove the original BAM file
            ;;
        "fastq")
            wget --continue --quiet -P "$output_dir" "$url"
            ;;
        "accession")
            prefetch "$url" --max-size 250G --output-directory "$output_dir"
            parallel-fastq-dump --sra-id "$url" --threads 32 --outdir "$output_dir" --split-files --gzip
            ;;
        *)
            echo "Unknown filetype: $filetype"
            ;;
    esac
}

export -f process_link

# Main script
main() {
    local data_dir="$1"
    local max_parallel_jobs="$2"

    # Create a temporary file to store all the jobs
    local job_file=$(mktemp)

    # Collect all jobs
    for obs_dir in "$data_dir"/*/; do
        local obs_id=$(basename "$obs_dir")
        local json_file="$obs_dir/observation.json"
        local fastq_dir="$obs_dir/fastq"

        if [[ -f "$json_file" ]]; then
            extract_links "$json_file" | while IFS=' ' read -r url filetype accession filename; do
                echo "process_link '$url' '$filetype' '$accession' '$filename' '$fastq_dir'" >> "$job_file"
            done
        else
            echo "observation.json not found in $obs_dir" >&2
        fi
    done

    # Run jobs in parallel
    parallel --jobs "$max_parallel_jobs" --line-buffer < "$job_file"

    # Clean up
    rm "$job_file"
}

# Check if required arguments are provided
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <data_directory> <max_parallel_jobs>"
    exit 1
fi

main "$1" "$2"
