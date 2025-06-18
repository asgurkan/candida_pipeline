import pandas as pd
import pytaxonkit

configfile: "config.yaml"

output = snakemake.output[0]
input_taxon = snakemake.input.blast
data_dir = snakemake.input.data_dir

blast_df = pd.read_csv(input_taxon, sep="\t", names=["qseqid","sseqid","pident","qcovs",
                                                     "length","bitscore","evalue","staxids",
                                                     "sskingdoms","sscinames","sblastnames"])

looking_taxas = blast_df["staxids"].unique()
result = pytaxonkit.lineage(looking_taxas, data_dir=data_dir)

result_short = result[["TaxID", "Lineage"]]

tax_level = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
for tax_index in range(len(tax_level)):
    result_short[tax_level[tax_index]] = result_short["Lineage"].apply(lambda x: x.split(";")[tax_index])

result_short.rename(columns={"TaxID":"staxids"}, inplace=True)

input_df_join = pd.merge(blast_df, result_short, how="left", left_on="staxids", right_on="staxids")
input_df_join.to_csv(output, sep="\t", index=False)