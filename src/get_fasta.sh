#!/bin/bash

# paths
centenarian_genomes="" # path to centenarian genomes
ad_genomes="" # path to AD genomes
home="" # path to working directory

#files
files=$(ls "path_to_merged_output") # path to hits_filter output

# analysis --> get fasta files of viral elements
for f in $files
do
	HG=$(basename $f | cut -f 1 -d '_')
	HG_path=$(ls $ad_genomes/$HG/$HG.hifi.hifiasm.p_ctg.fa)
    output_name=$(basename $f | cut -f 1 -d '.')
	echo $HG_path $f
	sbatch get_fasta_parallel.sh $HG_path $f output3/$output_name.fasta2
done

