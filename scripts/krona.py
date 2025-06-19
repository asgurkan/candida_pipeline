import os
import pandas as pd
import subprocess

# Get inputs and outputs from Snakemake
taxonomy_file = snakemake.input.taxonomy
shared_file = snakemake.input.shared
krona_html_output = snakemake.output.krona_html_output
krona_text_output = snakemake.output.krona_text_output

# Extract sample name
sample = os.path.basename(taxonomy_file).split(".")[0]

# Define output directories
TEXT_OUTPUT_DIR = os.path.dirname(krona_text_output)
HTML_BASE_DIR = os.path.dirname(krona_html_output)

# Ensure necessary directories exist
os.makedirs(TEXT_OUTPUT_DIR, exist_ok=True)
os.makedirs(HTML_BASE_DIR, exist_ok=True)

# Debugging: Print paths
print(f"Processing sample: {sample}")
print(f"Taxonomy file: {taxonomy_file}")
print(f"Shared file: {shared_file}")
print(f"Krona input file: {krona_text_output}")
print(f"Krona HTML output: {krona_html_output}")

# Read taxonomy data
taxonomy_df = pd.read_csv(taxonomy_file, sep="\t", dtype={"OTU": str})
taxonomy_df.columns = ["OTU", "Size", "Taxonomy"]

# Split taxonomy into separate columns
taxonomy_split = taxonomy_df["Taxonomy"].str.rstrip(";").str.split(";", expand=True)
taxonomy_split.columns = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
taxonomy_df = pd.concat([taxonomy_df["OTU"], taxonomy_split], axis=1)

# Read shared file
shared_df = pd.read_csv(shared_file, sep="\t", dtype={"OTU": str})

# Extract OTU abundance for the sample
otu_counts = shared_df.iloc[0, 3:].to_dict()  # Skip first 3 columns
otu_counts_df = pd.DataFrame(list(otu_counts.items()), columns=["OTU", "Count"])
otu_counts_df["OTU"] = otu_counts_df["OTU"].astype(str)
taxonomy_df["OTU"] = taxonomy_df["OTU"].astype(str)

# Merge data
merged_df = otu_counts_df.merge(taxonomy_df, on="OTU", how="inner")

# Keep only required columns for Krona
krona_df = merged_df[["Count", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]]
krona_df = krona_df.fillna("Unclassified")

# Ensure the output file is a file (not mistakenly a directory)
if os.path.isdir(krona_text_output):
    print(f"Error: {krona_text_output} exists as a directory, removing it...")
    os.rmdir(krona_text_output)  # Remove empty directory

# Save Krona input file
krona_df.to_csv(krona_text_output, sep="\t", index=False, header=False)

# Generate Krona plot
subprocess.run(["ktImportText", krona_text_output, "-o", krona_html_output])

print(f"Krona HTML report generated: {krona_html_output}")
