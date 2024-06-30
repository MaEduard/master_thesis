#!/bin/bash

# fasta files
files=$(ls "path_to_merged_output") # path to hits_filter output

# submit 1 reverse BLASTx job per file
for f in $files
do
	file_name=$(basename $f | cut -f 1 -d '.')
	sbatch reverse_diamond_parallel.sh $f output4/$file_name.reverse2.txt
done
