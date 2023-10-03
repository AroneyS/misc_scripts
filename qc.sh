#!/bin/bash

CONDA_ENV=bbtools
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $CONDA_ENV

JOBS=4
CPUS=16

while getopts 'i:o:c:j:' flag; do
  case "${flag}" in
    i) INPUT_FILE_LIST="${OPTARG}" ;;
    o) OUTPUT_DIR="${OPTARG}" ;;
    c) CPUS="${OPTARG}" ;;
    j) JOBS="${OPTARG}" ;;
    *) echo "Usage: -i input_file_list -o output_dir [-c cpus -j jobs]"
       exit 1 ;;
  esac
done

mkdir -p $OUTPUT_DIR/logs

cat $INPUT_FILE_LIST | parallel -j $JOBS --plus --col-sep "\t" \
  bbduk.sh \
    in1={1} in2={2} \
    out1=$OUTPUT_DIR/{1/} out2=$OUTPUT_DIR/{2/} \
    ref=adapters,phix \
    ktrim=r k=23 mink=11 hdist=1 \
    qtrim=r trimq=10 \
    minlen=30 \
    stats=$OUTPUT_DIR/logs/{1/...}.stats statscolumns=3 \
    t=$CPUS \
    '&>' $OUTPUT_DIR/logs/{1/...}.log
