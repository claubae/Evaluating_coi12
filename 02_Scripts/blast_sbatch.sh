#!/bin/bash
#SBATCH --job-name=blast_taxblock
#SBATCH --output=logs/blast_12_%j.out
#SBATCH --error=logs/blast_12_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

#Set environment
input_file=$1
input_taxid=$2
database=$3
blocked_taxID=$4
accession=$5 
output_folder=$6
cores=$SLURM_CPUS_PER_TASK   

startTime=$(date +%s)

outfile="${output_folder}/output_${input_taxid}_${accession}_blocked_${blocked_taxID}.txt"

echo "This script will run blast on $input_file (taxID: $input_taxid) with $cores cores, blocking this taxID: $blocked_taxID and write the output in $output_folder"

cd "$output_folder"


/home/cbaeyens/.conda/envs/blast/bin/blastn
export BLASTDB="/home/meg/dbs"
echo $(date +%H:%M) "BLASTing..."
cd $output_folder
blastn \
-query "$input_file" \
-db "$database" \
-num_threads "$cores" \
-perc_identity 80 \
-negative_taxids "$blocked_taxID" \
-word_size 11 \
-max_target_seqs 500 \
-outfmt "6 qseqid sseqid qlen slen pident mismatch evalue bitscore length staxids" \
-out "$outfile"

echo "BLAST finished. Result in: $outfile"
echo "It took $(expr $(date +%s) - $startTime) seconds"

echo "output file has "
wc -l "$outfile"
echo "lines"
echo
echo "It took $(expr `date +%s` - $startTime) seconds"
