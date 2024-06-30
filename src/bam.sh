#!/bin/bash

# paths
centenarian_genomes="" # path to centenarian genomes
ad_genomes="" # path to AD genomes

#files
files=$(ls *)

echo "Mapping files to hg38 genome..."

for f in $files
do
	g=$(basename $f | cut -f 1 -d '.' | cut -f 1 -d '_')
	file=$(basename $f | cut -f 1 -d '.')
	bam_file=$(ls $centenarian_genomes/$g/$g.hifi.hifiasm.p_ctg_hg38.bam)
	if [ -z "$bam_file" ]
	then
		echo "Cannot find the corresponding bam file"
	else
		echo "Mapping $file"
		sbatch bam_parallel.sh $bam_file $f output6/$file.mapped2.txt
	fi
done
