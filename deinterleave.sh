#!/bin/bash

PIGZ_COMPRESSION_THREADS=5

while getopts 'i:o:' flag; do
  case "${flag}" in
    i) INPUT_FILE_LIST="${OPTARG}" ;;
    o) OUTPUT_DIR="${OPTARG}" ;;
    *) echo "Usage: -i input_file_list -o output_dir"
       exit 1 ;;
  esac
done

mkdir -p $OUTPUT_DIR

cat $INPUT_FILE_LIST | parallel --plus \
    gunzip -c {} '|' \
    paste - - - - - - - - \
    '|' tee '>('cut -f 1-4 '|' tr '"\t"' '"\n"' '|' pigz --best --processes ${PIGZ_COMPRESSION_THREADS} '>' $OUTPUT_DIR/{/..}.1.fq.gz')' \
    '|' cut -f 5-8 '|' tr '"\t"' '"\n"' '|' pigz --best --processes ${PIGZ_COMPRESSION_THREADS} '>' $OUTPUT_DIR/{/..}.2.fq.gz > test.txt

