import pandas as pd 
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)

input = snakemake.input
qiime_output = snakemake.output.qiime_output
mothur_output = snakemake.output.mothur_output


df_taxon = pd.read_csv(str(input), sep="\t")

taxon_concat = lambda x:f"k__{x.Superkingdom}; p__{x.Phylum}; c__{x['Class']}; o__{x.Order}; f__{x.Family}; g__{x.Genus}; s__{x.Species}"
df_taxon["taxon"] = df_taxon.apply(lambda row: taxon_concat(row), axis=1)
df_taxon.rename(columns={"tax_id":"staxids"}, inplace=True)
df_taxon["estimated_counts"] = df_taxon["estimated_counts"].astype(int)

new_df = pd.DataFrame(columns=df_taxon.columns)  # Yeni bir boş DataFrame oluştur

# Her bir satırın çoğaltılması
for _, row in df_taxon.iterrows():
    tax_id = row['staxids']
    count = row['estimated_counts']
    for _ in range(count):
        new_df = new_df._append(row, ignore_index=True)


new_df = new_df[new_df["staxids"] != "unassigned"]
df_taxon_return = new_df[["staxids","taxon"]]
df_taxon_return.dropna(inplace=True)
df_taxon_return.to_csv(str(qiime_output), sep="\t",index=False)


#mothur_file_creating
df_mothur = pd.read_csv(str(input), sep="\t")
taxon_concat_mothur = lambda x:f"{x.Superkingdom}; {x.Phylum}; {x['Class']}; {x.Order}; {x.Family}; {x.Genus}; {x.Species}"
df_mothur["taxon"] = df_mothur.apply(lambda row: taxon_concat_mothur(row), axis=1)
df_mothur.rename(columns={"tax_id":"staxids"}, inplace=True)
df_mothur["estimated_counts"] = df_mothur["estimated_counts"].astype(int)

new_df2 = pd.DataFrame(columns=df_mothur.columns)  # Yeni bir boş DataFrame oluştur

# Her bir satırın çoğaltılması
for _, row in df_mothur.iterrows():
    tax_id = row['staxids']
    count = row['estimated_counts']
    for _ in range(count):
        new_df2 = new_df2._append(row, ignore_index=True)

new_df2 = new_df2[new_df2["staxids"] != "unassigned"]
df_mothur_return = new_df2[["staxids","taxon"]]
df_mothur_return.dropna(inplace=True)
df_mothur_return.to_csv(str(mothur_output), sep="\t",index=False)