#!/bin/bash
INPUT_FASTA=$1
OUTPUT_FASTA=$2
mafft --retree 2 --maxiterate 100 "$INPUT_FASTA"> "$OUTPUT_FASTA"