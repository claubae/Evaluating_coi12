#!usr/bin/bash

# First acitvate the conda environment with cutadapt

# Looks for the Leray Fwd sequence in a database, then removes it so all sequences are in the same frame
# usage bash <path_to_script> <fasta_file_db> <output_folder>
# Warning: The output folder must exist before executing the script.

LERAY_FWD=GGWACWGGWTGAACWGTWTAYCCYCC

INPUT_FASTA=$1

OUTPUT_FOLDER=$2

cutadapt -g ""${LERAY_FWD}";min_overlap=24" \
"${INPUT_FASTA}" -o "${OUTPUT_FOLDER}"/trimmed.fasta -e 4 -rc \
--discard-untrimmed \
--minimum-length 200