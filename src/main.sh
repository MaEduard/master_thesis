#!/bin/bash

# paths
centenarian_genomes="" # path to centenarian genomes
ad_genomes="" # path to AD genomes
home="" # path to working directory

# genomes
files=$(cat "my_path") # folder to genomes that need to be processed 

# databases
dbs=$(find databases/ -name "refseq_virus_unique.dmnd") # path to virus protein database

# date = id
today=$(date)
today=${today// /_}
today=${today//:/}

# start pipeline

for f in $files
do
	genomes=$(ls "path_to_fasta") # path to genome fasta file 
	for g in $genomes
	do
		for db in $dbs
		do
			seq_basename=$(basename $f)
			db_basename=$(basename $db)
			echo "$seq_basename"
			echo "$db_basename"
			echo $g
			sbatch diamond_forward_parallel.sh $db $g $today
		done
	done
done
