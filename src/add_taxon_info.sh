#!/bin/bash

files=$(find "path_to_BAM_analysis_output")

for f in $files
do
	echo $f
	dir_name=$(dirname $f)
	file_name=$(basename $f | cut -d '.' -f 1)
        python add_taxon_info.py $f refseq_unique_fasta_taxids.txt $dir_name/$file_name.fully_mapped2.txt
done
