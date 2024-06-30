import pytaxonkit
import pandas as pd
import sys

mapped_file = sys.argv[1]
acc2taxid_file = sys.argv[2]
outfile = sys.argv[3]

# Load data in dataframe
df = pd.read_csv(mapped_file, sep='\t')
acc2taxid = pd.read_csv(acc2taxid_file, names=["virus", "taxid"], sep='\t')

# Convert acc2taxid to dictionary
acc2taxid_dict =  dict(zip(acc2taxid.virus, acc2taxid.taxid))

# Create taxid column
df["taxid"] = df["accession_number"].str.replace(r'\.\d+', '', regex=True).map(acc2taxid_dict)

# Query taxonomy data from local database (obtained at ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz)
lineages = pytaxonkit.lineage(df["taxid"].to_list())

# Add taxonomy data to hit file
# df["superkingdom"], df["phylum"], df["class"], df["order"], df["family"], df["genus"], df["species"] = lineages["Lineage"].str.split(";").replace(r'^\s*$', 'X', regex=True)
lineages_split = lineages["Lineage"].str.split(";", expand=True)
lineages_split = lineages_split.fillna('X')  # Replace NaN values with 'X'
df["superkingdom"] = lineages_split[0]
df["phylum"] = lineages_split[1]
df["class"] = lineages_split[2]
df["order"] = lineages_split[3]
df["family"] = lineages_split[4]
df["genus"] = lineages_split[5]
df["species"] = lineages_split[6]

# Save updated hit file
df.to_csv(outfile, sep="\t", index=False)
