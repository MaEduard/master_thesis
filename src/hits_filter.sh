#!/bin/bash

# paths

# forward files
files=$(ls "path_to_forward_blastx_output") # path to the output of the foward BLASTx experiment

# start pipeline
for f in $files
do
	name=$(basename $f | cut -f 1 -d '.')
        sbatch hits_filter_parallel.sh $f output2/${name}.merged2.txt
done
