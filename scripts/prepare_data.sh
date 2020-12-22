#!/usr/bin/env bash

# usage: prepare_data.sh data_dir output.fasta quality_cutoff
# requirements: cutadapt

# check arguments
if [ "$#" -ne 2 -a "$#" -ne 3 ]; then
    echo "Error: wrong number of parameters"
    echo "Usage: prepare_data.sh data_dir output.fasta [quality_cutoff=30]"
    exit;
elif [ ! -d "$1" ]; then
    echo "Error: data directory $1 does not exist"
    exit;
fi
if [ -z "$3" ]; then
    echo "Quality cut-off not specified"
    echo "Using default quality cut-off value: 30"
fi

# input arguments
data_dir=$(realpath "$1")
out_fasta=$(realpath "$2")
q_cutoff=${3:-30}

outdir=$(dirname "$out_fasta")
out_filename=$(basename "$out_fasta" .fasta)
out_fastq="${outdir}/${out_filename}.fastq"

echo "----------------------------------------"
echo -e "input ab1 directory:\t" $data_dir
echo -e "output fasta:\t\t" $out_fasta
echo -e "quality cut-off:\t" $q_cutoff
echo "----------------------------------------"

# convert ab1 files to fastq
batch_ab1tofastq.py "$data_dir" "$out_fastq"

# trim sequences based on specified quality
echo "- trimming sequences ..."
cutadapt --quality-cutoff $q_cutoff,$q_cutoff --trim-n -m 100 --quality-base=64 \
    -o "$out_fasta" "$out_fastq" > "${outdir}/cutadapt.log"

echo "Done"
