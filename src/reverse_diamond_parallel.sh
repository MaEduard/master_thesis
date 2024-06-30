#!/bin/bash
#SBATCH -c 10

viral_insertions=$1
outputfile=$2
diamond blastx -d human_protein_database_nr/GRCh38_latest_protein_nr_db.dmnd -e 1e-06 --threads 9 -f 6 qseqid qstart qend salltitles evalue qframe pident qcovhsp sstart send slen -q $viral_insertions -o $outputfile
