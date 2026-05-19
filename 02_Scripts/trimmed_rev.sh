#!usr/bin/env bash

# First acitvate the conda environment with cutadapt (mirar fichita)
#Usage: bash Scripts/trimmed_rev_calanoida.sh input_file.fasta output_folder
#Looks for the Leray Rev sequence in a database, then removes it but keep the sequences that do not have it.
LERAY_RV=TGATTTTTTGGTCACCCTGAAGTTTA
INPUT_FASTA=$1
OUTPUT_FOLDER=$2

cutadapt \
  -a "${LERAY_RV};min_overlap=24" \
  -e 4 \
  -o "${OUTPUT_FOLDER}/trimmed_reverse.fasta" \
  "$INPUT_FASTA"
