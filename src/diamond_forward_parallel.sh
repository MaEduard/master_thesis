#!/bin/bash
#SBATCH -c 10

db=$1
db_name=${db%.dmnd}
db_name=$(basename $db_name)

HG=$2

genome_name_processed=${HG%.fa}
g_name=${genome_name_processed##*/}
g=${g_name//./_}

today=$3

time diamond blastx -d $db -e 1e-04 --threads 9 -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore salltitles qframe qcovhsp --sensitive -q $HG -o output/baseline_genomes_refseq/${g}-${db_name}_${today}.txt
