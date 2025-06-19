import os
import pandas as pd

blast_file = snakemake.input.blast_result
dbs_path = snakemake.params.dbs_folder
combined_result_text = snakemake.output.combined_result

fasta_folder = f'{dbs_path}/ref_genomes_clades/fna'

clade_ref_path_dict = {
    "clade1": f"{fasta_folder}/clade1_GCA_002759435.3.fna",
    "clade2": f"{fasta_folder}/clade2_GCF_003013715.1.fna",
    "clade3": f"{fasta_folder}/clade3_GCF_002775015.1.fna",
    "clade4": f"{fasta_folder}/clade4_GCA_003014415.1.fna",
    "clade5": f"{fasta_folder}/clade5_GCA_016809505.1.fna"
}

def get_fasta_headers_from_folder(fasta_folder):
    headers_dict = {}
    
    for fasta_file in os.listdir(fasta_folder):
        fasta_path = os.path.join(fasta_folder, fasta_file)

        if fasta_file.endswith('.fna'):
            with open(fasta_path, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        current_header = line.strip().split()[0][1:] 
                        clade = fasta_file.split('_')[0]  
                        headers_dict[current_header] = clade  
    return headers_dict

headers = get_fasta_headers_from_folder(fasta_folder)

def process_blast_results(blast_file, headers_dict):
    results_with_clade = []
    with open(blast_file, 'r') as f:
        for line in f:
            columns = line.strip().split('\t')  
            contig_id = columns[0]
            subject_id = columns[1]  

            if subject_id in headers_dict:
                clade_file = headers_dict[subject_id] 
                results_with_clade.append(columns + [clade_file])  
    return results_with_clade

results_with_clade = process_blast_results(blast_file, headers)

updated_blast_results = pd.DataFrame(results_with_clade, columns=[
    'Query ID', 'Subject ID', '% Identity', 'Alignment Length', 
    'Mismatches', 'Gap Openings', 'Q. Start', 'Q. End', 
    'S. Start', 'S. End', 'E-value', 'Bit Score', 'Clade'
])

updated_blast_results.iloc[:, 2:-1] = updated_blast_results.iloc[:, 2:-1].apply(pd.to_numeric, errors='coerce')

subset_50 = updated_blast_results[(updated_blast_results["E-value"] == 0) & (updated_blast_results["% Identity"] >= 90)] \
    .sort_values(by="Bit Score", ascending=False).head(50)

clade_counts = subset_50["Clade"].value_counts()
clade_percentage = (clade_counts / clade_counts.sum()) * 100

clade_percentage_with_paths = "\n".join([
    f"{clade} {percentage:.2f} {clade_ref_path_dict.get(clade, 'Path not found')}"
    for clade, percentage in clade_percentage.items()
])

with open(combined_result_text, "w") as f:
    f.write(clade_percentage_with_paths)