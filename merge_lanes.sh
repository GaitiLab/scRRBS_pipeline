#!/bin/bash

# The first argument is the output file path for the merged files.
OUTPUT_PATH=$1

# Remove the first argument from the list (it's the output path, not an input file)
shift

# Function to perform the merge operation
merge_files() {
    local output_path=$1
    shift  # Remove the output path argument

    # The rest of the arguments are input files to be merged
    local -a input_files=("$@")

    # Sort the input files
    IFS=$'\n' sorted_input_files=($(sort <<<"${input_files[*]}"))
    unset IFS

    # Determine the read type from the output file name
    if [[ "$output_path" =~ R1 ]]; then
        read_type="R1"
    elif [[ "$output_path" =~ R2 ]]; then
        read_type="R2"
    else
        echo "Error: Cannot determine read type from output file path."
        exit 1
    fi

    # Merge the input files
    echo "Merging $read_type"
    zcat "${sorted_input_files[@]}" > "${output_path%.fastq}.fastq"
    echo "Output file: ${output_path}"

    # Check if the output file exists
    if [[ ! -f "${output_path}" ]]; then
        echo "Error: Output file ${output_path} does not exist!"
        exit 1
    fi
}

# Call the merge_files function with the output path and input files
merge_files "$OUTPUT_PATH" "$@"
