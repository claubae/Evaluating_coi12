#Aquí partimos del supuesto que el marco de lectura correcto es 2, por tanto, cambiamos el valor de i (a i=2 para que elimine el 1er caracter, que corresponde a una tercera base), del módulo y también del archivo empleado:
#Usage: bash Scripts/Generate_coi12.sh input.fasta output.fasta
INPUT_FILE="$1"
OUTPUT_FILE="$2"
awk '/^>/ {print; next} {for (i=2; i<=length($0); i++) if (i % 3 != 1) printf "%s", substr($0, i, 1); printf "\n"}' "$INPUT_FILE" > "$OUTPUT_FILE"
#Primero se le dice que no haga nada (print;next) si la línea empieza por > (es decir, el encabezado)
#Luego que si la longitud de la secuencia no es 0, que vaya caracter por caracter calculando el módulo de su posición cuando se divide por 3-> si módulo no es igual a 1 (!0) se imprime (no se hace nada)
#Si modulo es igual a 1 simplemente ignora el caracter y NO lo imprime.