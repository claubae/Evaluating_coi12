#!/bin/bash
#SBATCH --job-name=blast_chunks
#SBATCH --output=logs/blast_12_%A_%a.out
#SBATCH --error=logs/error_12_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=08:00:00
#SBATCH --array=0-999             

/home/cbaeyens/.conda/envs/blast/bin/blastn 
export BLASTDB="/home/meg/dbs"
INPUT_DIR=$1
OUTPUT_DIR=$2
DATABASE="/home/meg/dbs/midori_12" 
FILES=(${INPUT_DIR}/*.fasta)   
TOTAL_FILES=${#FILES[@]}
NUM_TASKS=1000                    


FILES_PER_TASK=$(( (TOTAL_FILES + NUM_TASKS - 1) / NUM_TASKS ))
START_INDEX=$(( SLURM_ARRAY_TASK_ID * FILES_PER_TASK ))
END_INDEX=$(( START_INDEX + FILES_PER_TASK - 1 ))


if [ $END_INDEX -ge $TOTAL_FILES ]; then
    END_INDEX=$(( TOTAL_FILES - 1 ))
fi


for (( i=$START_INDEX; i<=$END_INDEX; i++ )); do
    QUERY_FILE="${FILES[$i]}"
    BASE_NAME=$(basename "$QUERY_FILE" .fasta)
    
    echo "Task $SLURM_ARRAY_TASK_ID processing file $i: $QUERY_FILE"

blastn \
      -query "$QUERY_FILE" \
      -db "$DATABASE" \
      -num_threads "$SLURM_CPUS_PER_TASK" \
      -perc_identity 90 \
      -word_size 11 \
      -max_target_seqs 500 \
      -outfmt "6 qseqid sseqid qlen slen pident mismatch evalue bitscore length staxids" \
      -out "${OUTPUT_DIR}/${BASE_NAME}_results.txt"
done
