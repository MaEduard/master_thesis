import pybedtools
import pandas as pd
import numpy as np
import sys

file = sys.argv[1]
outfile = sys.argv[2]

print(f"\n------ FILTERING {file} ------\n")

new_headers = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send slen evalue bitscore salltitles qframe qcovhsp"

df_hits = pd.read_csv(file, sep='\t', lineterminator='\n', names=new_headers.split())

print("done reading!")

if len(df_hits.index) == 0:
    print(f"0 hits found in {file}, stop processing")
else:
    df_hits.insert(len(df_hits.columns), "strandedness", ["" for i in range(len(df_hits.index))])

    file = file.split("/")
    file = file[len(file) - 1]
    file_splitted = file.split("_")
    virus_idx = ""
    if "protein" in file:
        virus_idx = file_splitted.index("protein") - 1
    else:
        virus_idx = file_splitted.index("virus") - 1
    virus_name = file_splitted[virus_idx]
    HG = file_splitted[0]

    df_hits["strandedness"] = df_hits.apply(lambda row: '-' if row["qstart"] > row["qend"] else '+', axis=1)

    mask = df_hits["qstart"] > df_hits["qend"]
    df_hits.loc[mask, ["qstart", "qend"]] = df_hits.loc[mask, ["qend", "qstart"]].values

    df_hits = df_hits.sort_values(["qseqid", "qstart", "qend"], ignore_index=True) # sort in order to merge
    df_hits.insert(len(df_hits.columns), "index", list(df_hits.index)) # keep track of what is merged through indices

    # reoder the columns to bed format
    df_hits = df_hits[["qseqid", "qstart", "qend", "sseqid", "bitscore", "strandedness", "pident", "length", "mismatch", "gapopen", "sstart", "send", "slen", "salltitles", "qframe", "qcovhsp","evalue", "index"]]
    x = pybedtools.BedTool.from_dataframe(df_hits)
    # merge hits based on overlapping genomic region
    c = x.merge(c=[17,18], o=["collapse","collapse"], s=True) # index in c is 1-based

    # only keep the hit in an overlapping region that has the highest bit score
    df_filtered_genomic_position_merged = []

    cnt = 0
    for i, line in enumerate(c):
        if cnt % 1000 == 0:
            print(f"Processing... {cnt/len(c) * 100}%")
        indices = line[4].split(',')
        indices = [int(float(i)) for i in indices]
        evalue = line[3].split(',')
        evalue = [float(i) for i in evalue]
        index_max_evalue = evalue.index(min(evalue))
        index_max_df = indices[index_max_evalue]
        df_filtered_genomic_position_merged.append(df_hits.iloc[index_max_df].tolist())
        cnt += 1

    df_filtered_genomic_position_merged = pd.DataFrame(df_filtered_genomic_position_merged, columns=["qseqid", "qstart", "qend", "sseqid", "bitscore", "strandedness", "pident", "length", "mismatch", "gapopen", "sstart", "send", "slen", "salltitles", "qframe", "qcovhsp", "evalue", "index"])
    print(f"Number of hits for {virus_name} in {HG}: {len(df_filtered_genomic_position_merged.index)} after merging overlapping hits with regions in the genome.")
    df_filtered_genomic_position_merged.to_csv(outfile, sep ='\t', header=True, index=False)

print(f"\n------ DONE FILTERING {file} ------\n")

