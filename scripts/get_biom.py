from biom.table import Table
from biom.util import biom_open

import pandas as pd


input_file = snakemake.input
sample_names = snakemake.params
print(sample_names)
output_file = snakemake.output

df = pd.read_csv(str(input_file), sep="\t")
sample_ids = sample_names

df_grouped = df.groupby("taxon").count().reset_index()
data = df_grouped["staxids"]

observ_ids = []
observ_metadata = []
for a_index, a_value in df_grouped.iterrows():
    tax_id = df[df["taxon"]==a_value["taxon"]]["staxids"].values[0]
    print(tax_id)
    observ_ids.append(str(tax_id))
    observ_metadata.append({'taxonomy': a_value["taxon"].split("; ")})

data_np = data.to_numpy()
data_np = data_np.reshape((len(data_np),1))

table = Table(data_np, observ_ids, sample_ids, observ_metadata, type="OTU table")

with biom_open(str(output_file), 'w') as f:
    table.to_hdf5(f, sample_names)