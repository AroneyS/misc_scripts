#!/bin/bash

while getopts 'i:o:n:' flag; do
  case "${flag}" in
    i) INPUT_DIR="${OPTARG}" ;;
    o) OUTPUT_DIR="${OPTARG}" ;;
    n) BATCH_SIZE="${OPTARG}" ;;
    *) echo "Usage: -i input_dir -o output_dir -n batch_size"
       exit 1 ;;
  esac
done

find $INPUT_DIR -maxdepth 1 -type f |
  parallel -n$BATCH_SIZE \
  mkdir $OUTPUT_DIR/batch{#} \
  ';' ln -s {} $OUTPUT_DIR/batch{#}
