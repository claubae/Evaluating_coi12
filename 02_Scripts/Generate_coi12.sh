#Aquí partimos del supuesto que el marco de lectura correcto es 2, por tanto, cambiamos el valor de i, del módulo y también del archivo empleado:
#Usage: bash Scripts/Generate_coi12.sh input.fasta output.fasta
INPUT_FILE="$1"
OUTPUT_FILE="$2"
awk '/^>/ {print; next} {for (i=2; i<=length($0); i++) if (i % 3 != 1) printf "%s", substr($0, i, 1); printf "\n"}' "$INPUT_FILE" > "$OUTPUT_FILE"