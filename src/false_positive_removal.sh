#!/bin/bash

# extract hit files
files=$(ls "path_to_merged_output") # path to hits_filter output

# remove false positives
for f in $files
do
	file_name=$(basename $f | cut -f 1 -d '.')
	reverse_file=$(ls output4/$file_name.reverse2.filtered.txt)
	echo $reverse_file
	sbatch false_positive_removal_parallel.sh $f $reverse_file output5/$file_name
done
