# imports
import sys
import pandas as pd
import pybedtools

# files
hit_list_file = sys.argv[1]
reverse_diamond_file = sys.argv[2]
output_file = sys.argv[3]

headers_reverse_diamond = "qseqid qstart qend salltitles evalue qframe pident qcovhsp sstart send slen"
reverse_diamond_df = pd.read_csv(reverse_diamond_file, sep="\t", names=headers_reverse_diamond.split())
reverse_diamond_df.insert(len(reverse_diamond_df.columns), "strandedness", ["" for i in range(len(reverse_diamond_df.index))])

if len(reverse_diamond_df) == 0:
    content = ""
    with open(hit_list_file, 'r') as source_file:
        content = source_file.read()

    with open(output_file, 'w') as file:
        file.write(content)
        print("No false positives found")
else:
    # Add strandedness column
    reverse_diamond_df["strandedness"] = reverse_diamond_df.apply(lambda row: '-' if row["qstart"] > row["qend"] else '+', axis=1)

    # Switch the order of start and stop positions for entries with reversed hits
    reverse_entries = reverse_diamond_df["qstart"] > reverse_diamond_df["qend"]
    reverse_diamond_df.loc[reverse_entries, ["qstart", "qend"]] = reverse_diamond_df.loc[reverse_entries, ["qend", "qstart"]].values

    # Split the 'qseqid' column into two parts based on ':'
    reverse_diamond_df[['seqid_part', 'coordinates']] = reverse_diamond_df['qseqid'].str.split(':', expand=True)

    # Split the 'coordinates' column into two parts based on '-'
    reverse_diamond_df[['start_coord', 'end_coord']] = reverse_diamond_df['coordinates'].str.split('-', expand=True)

    # Convert the new columns to numeric values
    reverse_diamond_df['start_coord'] = pd.to_numeric(reverse_diamond_df['start_coord'])
    reverse_diamond_df['end_coord'] = pd.to_numeric(reverse_diamond_df['end_coord'])

    # Update 'qstart' and 'qend' columns
    reverse_diamond_df['qstart'] = reverse_diamond_df['qstart'] + reverse_diamond_df['start_coord']
    reverse_diamond_df['qend'] = reverse_diamond_df['qend'] + reverse_diamond_df['end_coord']

    # Drop the intermediate columns if needed
    reverse_diamond_df.drop(['qseqid', 'coordinates', 'start_coord', 'end_coord'], axis=1, inplace=True)

    # Convert to BedTool
    reverse_diamond_bedtool = pybedtools.BedTool.from_dataframe(reverse_diamond_df["seqid_part qstart qend salltitles evalue strandedness".split()])

    hit_list_df = pd.read_csv(hit_list_file, sep="\t")

    hit_list_bedtool_intersect =  pybedtools.BedTool.from_dataframe(hit_list_df)

    hit_list_df["old_qstart"] = hit_list_df["qstart"]
    hit_list_df["old_qend"] = hit_list_df["qend"]

    hit_list_bedtool_substract =  pybedtools.BedTool.from_dataframe(hit_list_df)

    d = hit_list_bedtool_intersect.intersect(reverse_diamond_bedtool, s=True, v=True, f=0.5)

    d_name = output_file + ".intersect2.txt"
    d.moveto(d_name)

    print(f"Done processing {hit_list_file}!")
