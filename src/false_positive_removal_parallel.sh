#!/bin/bash
#SBATCH -c 1

hit_file=$1
reverse_diamond_file=$2
output_file=$3

echo "Processing: $1 $2 $3"
python false_positive_removal.py $1 $2 $3
