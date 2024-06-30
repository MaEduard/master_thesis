#!/bin/bash
#SBATCH -c 1

bam_file=$1
f=$2
output=$3

python bam_analysis_updated.py $bam_file $f $output
