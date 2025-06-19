import pandas as pd
from shutil import copytree
import my_plots as mp
from krona_plots import krona_report
#from create_word import create_word_file, move_fastqcs
from create_html_report import create_html_report
import os

def get_graph_path(path_str):
    path_str_split = path_str.split("/")
    return "/".join(path_str_split[3:])


output = snakemake.output.excel
project_name = output.split("/")[-2]
input_summary = snakemake.input.summary
metadata = snakemake.input.metadata
beta = snakemake.input.beta
alpha = snakemake.input.alpha
css_style = snakemake.params.css_style
organism = snakemake.params.organism


output_files_dict = dict()



# input_summary = "data/results/ind_results/merged.tax.summary"
# metadata = "medical_metadata.tsv"
# beta="data/results/ind_results/beta-diversity/beta_diversity_data.xlsx"
# alpha="data/results/ind_results/alpha-diversity/alpha_diversity_data.tsv"

# output = "data/report/Excel_Report.xlsx"


summary=pd.read_csv(input_summary, sep="\t")
metadata_df = pd.read_csv(metadata, sep="\t")
metadata_df = metadata_df[metadata_df["proje"]==project_name]
output_files_dict["metadata"]=metadata_df

metadata_df_sample_dict = {b_value["sample-id"]:b_value["sample"] for a_key,b_value in metadata_df.iterrows()}

for csc_name in summary.columns:
    if csc_name in metadata_df_sample_dict.keys():
        summary.rename(columns={csc_name:metadata_df_sample_dict[csc_name]},inplace=True)


summary=summary[summary["taxon"]!="Root"]
tax_levels={1:"Domain", 2:"Phylum", 3:"Class", 4:"Order", 5:"Family", 6:"Genus", 7:"Species"}
tax_levels_base={0:"Zero", 1:"Domain", 2:"Phylum", 3:"Class", 4:"Order", 5:"Family", 6:"Genus", 7:"Species"}
summary["Taxonomy"] = summary["taxlevel"].map(lambda x: tax_levels[x])
summary["Mean"] = summary["total"]
summary.drop(columns=["daughterlevels","taxlevel","total"], inplace=True)
summary.reset_index(inplace=True)
summary["taxon"] = summary["taxon"].str.rstrip()

metadata_df = pd.read_csv(metadata, sep="\t")
metadata_df = metadata_df[metadata_df["proje"]==project_name]
len_stability = metadata_df.shape[0]
summary["Mean"] = summary["Mean"].apply(lambda x: x/len_stability)
summary["Mean"] = summary["Mean"].round(3)

alpha_df = pd.read_csv(alpha, sep="\t", index_col=0).rename_axis('sample-id').reset_index()
alpha_j_df = alpha_df.merge(metadata_df, on="sample-id", how="left")
alpha_j_df.sort_values(by="sample", inplace=True)

writer = pd.ExcelFile(str(beta))

jaccard_df = pd.read_excel(writer, sheet_name='jaccard', index_col=0).rename_axis('sample-id').reset_index()
jaccard_meta_df = jaccard_df.merge(metadata_df, on="sample-id", how="left")

bray_df = pd.read_excel(writer, sheet_name='bray_curtis', index_col=0).rename_axis('sample-id').reset_index()
bray_meta_df = bray_df.merge(metadata_df, on="sample-id", how="left")

unweighted_unifrac_df = pd.read_excel(writer, sheet_name='unweighted_unifrac', index_col=0).rename_axis('sample-id').reset_index()
unweighted_unifrac_df = unweighted_unifrac_df.merge(metadata_df, on="sample-id", how="left")

weighted_unifrac_df = pd.read_excel(writer, sheet_name='weighted_unifrac', index_col=0).rename_axis('sample-id').reset_index()
weighted_unifrac_meta_df = weighted_unifrac_df.merge(metadata_df, on="sample-id", how="left")


for i_index,i_value in metadata_df.iterrows():
    summary.rename(columns={i_value["sample-id"]:i_value["sample"]}, inplace=True)

    alpha_j_df.rename(columns={i_value["sample-id"]:i_value["sample"]}, inplace=True)
    alpha_j_df["sample-id"] = alpha_j_df["sample-id"].replace([i_value["sample-id"]],i_value["sample"])

    jaccard_df.rename(columns={i_value["sample-id"]:i_value["sample"]}, inplace=True)
    jaccard_df["sample-id"] = jaccard_df["sample-id"].replace([i_value["sample-id"]],i_value["sample"])

    bray_df.rename(columns={i_value["sample-id"]:i_value["sample"]}, inplace=True)
    bray_df["sample-id"] = bray_df["sample-id"].replace([i_value["sample-id"]],i_value["sample"])

    unweighted_unifrac_df.rename(columns={i_value["sample-id"]:i_value["sample"]}, inplace=True)
    unweighted_unifrac_df["sample-id"] = unweighted_unifrac_df["sample-id"].replace([i_value["sample-id"]],i_value["sample"])

    weighted_unifrac_df.rename(columns={i_value["sample-id"]:i_value["sample"]}, inplace=True)
    weighted_unifrac_df["sample-id"] = weighted_unifrac_df["sample-id"].replace([i_value["sample-id"]],i_value["sample"])

def give_abundance(x, column, base_row_df):
    x_domain = x["rankID"][:3]
    base_row_df_ = base_row_df[base_row_df["rankID"] == x_domain]
    if base_row_df_[column].values[0] == 0:
        return 0

    abundance = x[column]*100/base_row_df_[column].values[0]

    return abundance.round(3)

def summary_func(summary_):
    base_columns = ["index","rankID","taxon","Taxonomy"]
    sample_columns = [a for a in summary_.columns if a not in base_columns]

    remove_columns=[]
    new_sample_columns=[]
    for sample_column in sample_columns:
        if sample_column!="Mean":
            remove_columns.append(sample_column)
            column_count = sample_column + "_Count"

            summary_[column_count] = summary[sample_column]
            new_sample_columns.append(sample_column)

    new_sample_columns.append("Mean")
    summary_["Mean_Count"] = summary_["Mean"]
    base_row = summary_[summary_["Taxonomy"]=="Domain"]
    for column in new_sample_columns:
        summary_[column] = summary_.apply(lambda x: give_abundance(x, column, base_row), axis=1)

    return_cols = base_columns + new_sample_columns

    summary_[new_sample_columns] = summary_[new_sample_columns].astype(float)
    return summary_, new_sample_columns

for s_index, s_item in summary.iterrows():
    if s_item["taxon"].startswith("unknown"):
        new_rank_id = ".".join(s_item["rankID"].split(".")[:-1])
        new_taxon_name = summary[summary["rankID"]==new_rank_id]["taxon"].values[0].rstrip()
        summary.loc[summary["rankID"]==s_item["rankID"], ["taxon"]] = f"{new_taxon_name} {s_item['Taxonomy']}"


summary_1 = summary.copy()

s_summary, sample_columns = summary_func(summary_1)
sample_columns.sort()

if organism == "Bacteria":
    bacteria_rankid = s_summary.query("taxon=='Bacteria'")["rankID"].values[0]
else: 
    bacteria_rankid = s_summary.query("taxon=='Eukaryota'")["rankID"].values[0]

summary_1_genus = s_summary.query("Taxonomy == 'Genus'")
summary_1.to_csv("summary_1.csv",sep=",")
krona_output_path =krona_report(summary_1, metadata_df, output) #krona report part
output_files_dict["krona"]=krona_output_path


writer.close()

# create_word_file(metadata_df,output,s_summary_sheet_phylum)
# move_fastqcs(output)
create_html_report(output_files_dict, output)
print(css_style)
copytree(css_style, f"data/report/{project_name}/styles")
