# imports
import glob
import pandas as pd
import pybedtools

#repeats = pd.read_csv("grch38.rmsk.bed", sep="\t|;", names=["ref", "ref_start", "ref_end", "program", "class", "family/motif"], engine="python")
repeats = pd.read_csv("data/rmsk_excl_ERVs.txt", sep="\t", header=None)

# Loading data
# Notice that the format of the file is very specific to the 'primary' aligment batches that were processed for refseq database.

df = pd.DataFrame({"contig": [],
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
                       "strandedness_ref": [],
                       "taxid": [],
                       "superkingdom": [],
                       "phylum": [],
                       "class": [],
                       "order": [],
                       "family": [],
                       "genus": [],
                       "species": []})

path_to_cohort_AD1 = "../mapped_files2/ad1_refseq/*fully_mapped2*"
path_to_cohort_AD2 = "../mapped_files2/ad2_refseq/*fully_mapped2*"
path_to_cohort_AD3 = "../mapped_files2/ad3_refseq/*fully_mapped2*"
path_to_cohort_AD4 = "../mapped_files2/ad4_refseq/*fully_mapped2*"
path_to_cohort_AD5 = "../mapped_files2/ad5_refseq/*fully_mapped2*"
path_to_cohort_cent1 = "../mapped_files2/centenarian1_refseq/*fully_mapped2*"
path_to_cohort_cent2 = "../mapped_files2/centenarian2_refseq/*fully_mapped2*"
path_to_cohort_cent3 = "../mapped_files2/centenarian3_refseq/*fully_mapped2*"
path_to_cohort_cent4 = "../mapped_files2/centenarian4_refseq/*fully_mapped2*"
path_to_cohort_cent5 = "../mapped_files2/centenarian5_refseq/*fully_mapped2*"
path_to_cohort_cent6 = "../mapped_files2/centenarian6_refseq/*fully_mapped2*"

# load files
files_AD1 = glob.glob(str(path_to_cohort_AD1))
files_AD2 = glob.glob(str(path_to_cohort_AD2))
files_AD3 = glob.glob(str(path_to_cohort_AD3))
files_AD4 = glob.glob(str(path_to_cohort_AD4))
files_AD5 = glob.glob(str(path_to_cohort_AD5))
files_cent1 = glob.glob(str(path_to_cohort_cent1))
files_cent2 = glob.glob(str(path_to_cohort_cent2))
files_cent3 = glob.glob(str(path_to_cohort_cent3))
files_cent4 = glob.glob(str(path_to_cohort_cent4))
files_cent5 = glob.glob(str(path_to_cohort_cent5))
files_cent6 = glob.glob(str(path_to_cohort_cent6))

all_files = files_AD1 + files_AD2 + files_AD3 + files_AD4 + files_AD5 + files_cent1 + files_cent2 + files_cent3 + files_cent4 + files_cent5 + files_cent6

AD_files = set(files_AD1 + files_AD2 + files_AD3 + files_AD4 + files_AD5)

for file in all_files:
    df_virus = pd.read_csv(file, sep="\t")
    sample = file.split('/')[-1].split('_')[0]

    cohort =  ""
    if file in AD_files:
        cohort = "AD"
    else:
        cohort = "Centenarian"

    # Add virus name column and sample column
    df_virus.insert(len(df_virus.columns), "sample", [sample]*len(df_virus.index), False)
    df_virus.insert(len(df_virus.columns), "cohort", [cohort]*len(df_virus.index), False)
    df = pd.concat([df, df_virus], ignore_index=True)

# Calculate insertion lenghts
df["length"] = pd.to_numeric(df["contig_end"]) - pd.to_numeric(df["contig_start"])

# Do not take unmapped integrations into account for now --> but might be interesting to look at still. They are just not mapped to the human genome but found in our contigs. 
df = df[~df["contig"].str.contains("_UNMAPPED")]

# See how much is left after removing all everything that slightly overlaps with the repeat masker
repeats = repeats[[repeats.columns[5], repeats.columns[6], repeats.columns[7], repeats.columns[10], repeats.columns[11], repeats.columns[9], repeats.columns[12]]]
repeat_bed = pybedtools.BedTool.from_dataframe(repeats)

df["ref_start"] = df["ref_start"].astype(int)
df["ref_end"] = df["ref_end"].astype(int)

df = df[["ref", "ref_start", "ref_end", "accession_number", "evalue", "strandedness_contig",
"virus_prot_start", "virus_prot_end", "contig", "contig_start", "contig_end", "pident", "mismatch",
"qframe", "contig_sequence", "cigar_so_far", "salltitles", "strandedness_ref", "taxid", "superkingdom",
"phylum", "class", "order", "family", "genus", "species", "cohort", "sample", "length"]]
df_bed = pybedtools.BedTool.from_dataframe(df)

print("Filtering!")
df_filtered = df_bed.intersect(repeat_bed, v=True, s=True, f=0.5)

df_filtered.moveto("final_dataset.txt")
