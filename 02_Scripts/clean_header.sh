#!usr/bin/env bash
#Usage: bash Scripts/clean_header input_file.fasta output_file.fasta
#Lo que quiero hacer con este script es limpiar el header: mantener solamente el accession number y el taxid
INPUT_FILE="$1"
OUTPUT_FILE="$2"
awk -F'[.;_]' '
    /^>/ {
        print $1, $NF
    } 
    !/^>/ {
        print $0
    }
' OFS=";taxid=" "$INPUT_FILE" > "$OUTPUT_FILE"