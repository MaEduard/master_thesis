#!/bin/bash
#SBATCH -c 1

echo "Hits filtering: $1 & $2"
python hits_filter.py $1 $2
