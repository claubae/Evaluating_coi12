#Este script es para generar una base de datos local para BLAST
#Usage: First activate the blast environment. bash Scripts/Makelocaldb.sh input.fasta output.file
INPUT_FILE="$1"
OUTPUT_FILE="$2" 
makeblastdb -in "$INPUT_FILE" \
            -dbtype nucl \
            -out "$OUTPUT_FILE" \
            -parse_seqids \
            -taxid_map "Prueba/taxmap.txt"