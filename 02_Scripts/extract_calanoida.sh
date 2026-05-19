#Usage: First activate de environment for blast. A escribir en la terminal: bash Scripts/extract_calanoida.sh path to the database (indicate the basename of the db) output.fasta.
#Make sure to download the NCBI taxonomy file first. Pon dichos archivos en la carpeta "raíz", no en la de la db (o cambiar de directorio antes de ejecutar el script). 
INPUT_FILE="$1"
OUTPUT_FILE="$2"
export BLASTDB="Downloads/TFG_COI/Prueba/dbs"
blastdbcmd -db "$INPUT_FILE" \
           -taxids 6833 \
           -outfmt %T,%s,%a,%S \
           -out "$OUTPUT_FILE" \