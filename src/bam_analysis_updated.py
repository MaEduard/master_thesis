import pysam
import pandas as pd
import sys
import re

# Maps the coordinate positions in the final tab-separated output file to positions in the human reference genome hg38.

# Args:
# bam_file_path: bam file containing contig -> hg38 mappings
# hit_file: .txt tab separated file containing processed human genome virus alignments
# output_file: .txt output file

bam_file_path = sys.argv[1]
hit_file = sys.argv[2]
numb_cols = 0
output_file = sys.argv[3]

output_df = pd.DataFrame({"contig": [],
                       "contig_start": [],
                       "contig_end": [],
                       "accession_number": [],
                       "virus_prot_start": [],
                       "virus_prot_end": [],
                       "strandedness_contig": [],
                       "pident": [],
                       "mismatch": [],
                       "qframe": [],
                       "ref": [],
                       "ref_start": [],
                       "ref_end": [],
                       "contig_sequence": [],
                       "cigar_so_far": [],
                       "salltitles": [],
                       "evalue": [],
                       "strandedness_ref": []})


# Convert bam AlignmentFile from pysam to dataframe including only columns of interest (i.e. alignment start and end in both contig and reference)

bam_df = pd.DataFrame({"contig": [],
                       "contig_start": [],
                       "contig_end": [],
                       "ref": [],
                       "ref_start": [],
                       "ref_end": [],
                       "CIGAR": [],
                       "query_sequence": [],
                       "is_supplementary":[],
                       "is_secondary": []})

def parse_cigar(cigar, start_virus_contig, end_virus_contig, start_ref):
    current_position_contig = 0
    current_position_ref = start_ref
    start_final = 0
    end_final = 0
    found_start = False
    operations = [op for op in re.split(r'(\D+)', cigar) if op]

    cigar_so_far = ""

    i = 0
    while i < len(operations):
        length = int(operations[i])
        op_type = operations[i+1]

        if op_type == 'M' or op_type == '=' or op_type == 'X':
            current_position_ref += length
            current_position_contig += length
        elif op_type == 'D' or op_type == 'N':
            current_position_ref += length
        elif op_type == 'I' or op_type == 'S' or op_type == 'H':
            current_position_contig += length
        else:
            raise ValueError(f"Symbol unknown: {op_type}")

        cigar_so_far += str(length) + op_type

        if not found_start and start_virus_contig <= current_position_contig:
            found_start = True

            if op_type == 'M' or op_type == '=' or op_type == 'X' or op_type == 'D' or op_type == 'N':
                start_final = current_position_ref - (current_position_contig - start_virus_contig)
            elif op_type == 'I' or op_type == 'S' or op_type == 'H':
                start_final = current_position_ref

        if found_start and end_virus_contig <= current_position_contig:

            if op_type == 'M' or op_type == '=' or op_type == 'X' or op_type == 'D' or op_type == 'N':
                end_final = current_position_ref - (current_position_contig - end_virus_contig)
            elif op_type == 'I' or op_type == 'S' or op_type == 'H':
                end_final = current_position_ref

            return start_final, end_final, cigar_so_far

        i += 2

    return None

line = ""
with open(hit_file, 'r') as temp:
    line = temp.readline()
    numb_cols = len(line.split('\t'))

headers = []

if numb_cols == 18:
    headers = ['qseqid', 'qstart', 'qend', 'sseqid', 'bitscore', 'strandedness', 'pident', 'length', 'mismatch', 'gapopen', 'sstart', 'send', 'slen', 'salltitles', 'qframe', 'qcovhsp', 'evalue', 'index']
elif numb_cols == 1:
    output_df.to_csv(output_file, sep="\t")
    raise Exception(f"Input {hit_file} is empty!")
else:
    print(line)
    print(numb_cols)
    raise Exception(f"Input {hit_file} not recognized based on the number of columns")

hits = pd.read_csv(hit_file, sep='\t', lineterminator='\n', names=headers)
unique_contig_names = set(hits["qseqid"])
to_save_alignmentSegments = {contig: [] for contig in unique_contig_names}

# Obtain all alignments of interest
print("Obtaining alignments of interest")
with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
    for alignment in bam_file:
        if alignment.query_name in unique_contig_names:
            to_save_alignmentSegments[alignment.query_name].append(alignment)

# Iterate over and process all EVEs one-by-one
print("Starting to map back the contig locations to hg38")
for i in range(len(hits.index)):
    query = hits.iloc[i]["qseqid"]
    is_reverse = hits.iloc[i]["strandedness"]
    alignments = to_save_alignmentSegments[query]

    if len(alignments) == 0:
        print(f"no alignment found for {query}")
    else:
        found_one = False
        for alignment in alignments:
            if not alignment.is_secondary and alignment.cigarstring != None:
                start_virus_HG, end_virus_HG, cigar_so_far = parse_cigar(alignment.cigarstring, hits.iloc[i]["qstart"], hits.iloc[i]["qend"],  alignment.reference_start)
                strandedness = ""
                if alignment.is_reverse:
                    strandedness = "+"
                else:
                    strandedness = "-"
                found_one = True
                output_df.loc[len(output_df.index)] = [hits.iloc[i]["qseqid"], hits.iloc[i]["qstart"], hits.iloc[i]["qend"], hits.iloc[i]["sseqid"], hits.iloc[i]["sstart"], hits.iloc[i]["send"], hits.iloc[i]["strandedness"], hits.iloc[i]["pident"], hits.iloc[i]["mismatch"], hits.iloc[i]["qframe"], alignment.reference_name , start_virus_HG, end_virus_HG, alignment.query_sequence[hits.iloc[i]["qstart"] - 1: hits.iloc[i]["qend"]], cigar_so_far, hits.iloc[i]["salltitles"], hits.iloc[i]["evalue"], strandedness]
                break
            
        if not found_one:
            mq = 0
            for alignment in alignments:
                if alignment.mapping_quality > mq and alignment.cigarstring != None:
                    mq = alignment.mapping_quality
                    start_virus_HG, end_virus_HG, cigar_so_far = parse_cigar(alignment.cigarstring, hits.iloc[i]["qstart"], hits.iloc[i]["qend"],  alignment.reference_start)
                    strandedness = ""
                    if alignment.is_reverse:
                        strandedness = "+"
                    else:
                        strandedness = "-"
                    output_df.loc[len(output_df.index)] = [hits.iloc[i]["qseqid"], hits.iloc[i]["qstart"], hits.iloc[i]["qend"], hits.iloc[i]["sseqid"], hits.iloc[i]["sstart"], hits.iloc[i]["send"], hits.iloc[i]["strandedness"], hits.iloc[i]["pident"], hits.iloc[i]["mismatch"], hits.iloc[i]["qframe"], alignment.reference_name , start_virus_HG, end_virus_HG, alignment.query_sequence[hits.iloc[i]["qstart"] - 1: hits.iloc[i]["qend"]], cigar_so_far, hits.iloc[i]["salltitles"], hits.iloc[i]["evalue"], strandedness]
                    found_one = True
            
            if not found_one:
                output_df.loc[len(output_df.index)] = [hits.iloc[i]["qseqid"] + "_UNMAPPED", hits.iloc[i]["qstart"], hits.iloc[i]["qend"], hits.iloc[i]["sseqid"], hits.iloc[i]["sstart"], hits.iloc[i]["send"], hits.iloc[i]["strandedness"], hits.iloc[i]["pident"], hits.iloc[i]["mismatch"], hits.iloc[i]["qframe"], 'X' , 'X', 'X', 'X', 'X', hits.iloc[i]["salltitles"], hits.iloc[i]["evalue"], 'X']                
                print(f"no alignment found for {query}")

print(f"done mappnig contig locations of {hit_file} to hg38")
output_df.to_csv(output_file, sep="\t", index=False)

