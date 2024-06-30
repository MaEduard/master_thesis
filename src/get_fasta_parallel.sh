#!/bin/bash
#SBATCH -c 1

echo "Processing get_fasta.sh with: $1 & $2 & $3"
name=$(basename $1)
file_name=$(basename $2)
tail -n +2 $2 > output2/merged_tmp_$file_name.txt
mkdir human_genomes/tmp_$file_name
cp $1 human_genomes/tmp_$file_name/

bedtools getfasta -fi human_genomes/tmp_$file_name/$name -bed output2/merged_tmp_$file_name.txt -fo $3

rm output2/merged_tmp_$file_name.txt
rm human_genomes/tmp_$file_name/$name*
rmdir human_genomes/tmp_$file_name

