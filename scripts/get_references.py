import pandas as pd

input_rel_abundance = snakemake.input.rel_abundance
database = snakemake.input.database
output = snakemake.output

input_summary_df = pd.read_csv(str(input_rel_abundance), sep="\t")
input_summary_df["tax_id"] = input_summary_df["tax_id"].astype(str)

def parse_fasta(file_path):
    headers = []
    sequences = []
    
    with open(file_path, "r") as file:
        header = ""
        sequence = ""
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                # Başlık satırı
                if header != "":
                    headers.append(header)
                    sequences.append(sequence)
                header = line[1:]
                sequence = ""
            else:
                # Sekans satırı
                sequence += line
        # Son sekansı ekle
        headers.append(header)
        sequences.append(sequence)
    
    # DataFrame oluşturma
    df = pd.DataFrame({"Header": headers, "Sequence": sequences})
    return df

fasta_df = parse_fasta(str(database))
fasta_df["Taxid"] = fasta_df["Header"].str.split(":", expand=True)[0]
fasta_df["Taxid"] = fasta_df["Taxid"].astype(str)

input_summary_df_with_sequence= pd.merge(input_summary_df, fasta_df, how="left", right_on="Taxid",left_on="tax_id")
input_summary_df_with_sequence = input_summary_df_with_sequence[["Taxid","Sequence"]]
input_summary_df_with_sequence.dropna(inplace=True)
input_summary_df_with_sequence['Taxid'] = ">" + input_summary_df_with_sequence["Taxid"].astype(str)

input_summary_df_with_sequence.to_csv(str(output), sep="\t", index=False, header=False)