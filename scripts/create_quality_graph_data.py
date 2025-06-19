import pandas as pd
from snakemake import shell

# inputs
BEFORE_FASTQC = snakemake.input.before_fastqc
AFTER_FASTQC = snakemake.input.after_fastqc

# outputs
BEFORE_PER_SEQ_QUALITY_SCORES = snakemake.output.before_per_seq_quality_scores
BEFORE_GENERAL_STATS = snakemake.output.before_general_stats

AFTER_PER_SEQ_QUALITY_SCORES = snakemake.output.after_per_seq_quality_scores
AFTER_GENERAL_STATS = snakemake.output.after_general_stats

# params
SAMPLE_NAME = snakemake.params.sample_name

# functions
def fetch_data(path, start_marker, end_marker):
    data = []
    with open(path, "r") as file:
        k = False
        for line in file:
            line = line.strip()
            if line.startswith(str(start_marker)):
                k = True
            if line.startswith(str(end_marker)):
                k = False
            if k == True:
                data.append(line)
    return data

def unzip(fastqc):
    zip_path = f"{fastqc}/{SAMPLE_NAME}_fastqc.zip"
    unzip_command = f"unzip {zip_path}  {SAMPLE_NAME}_fastqc/fastqc_data.txt -d {fastqc}"
    shell(unzip_command)

unzip(BEFORE_FASTQC)
unzip(AFTER_FASTQC)

def create_fastqc_raw_data_path(fastqc):
    fastqc_raw_data_path = f"{fastqc}/{SAMPLE_NAME}_fastqc/fastqc_data.txt"
    return fastqc_raw_data_path


def generate_general_stats_data(fastqc, general_stats_data_path):
    fastqc_raw_data_path = create_fastqc_raw_data_path(fastqc)
    general_stats_data = fetch_data(fastqc_raw_data_path, ">>Basic Statistics", ">>END_MODULE")
    general_stats_data = general_stats_data[2:]

    general_stats_data = pd.DataFrame([x.split("\t") for x in general_stats_data], columns=["Measure", "Value"])
    specific_measures = ["Total Sequences", "Total Bases", "Sequence length", "%GC"]
    general_stats_data = general_stats_data[general_stats_data["Measure"].isin(specific_measures)]
    general_stats_data.to_csv(general_stats_data_path, sep="\t", index=False)


def generate_per_seq_quality_scores_data(fastqc, per_seq_quality_scores_data_path):
    fastqc_raw_data_path = create_fastqc_raw_data_path(fastqc)

    per_seq_quality_scores_data = fetch_data(fastqc_raw_data_path, ">>Per sequence quality scores", ">>END_MODULE")
    per_seq_quality_scores_data = per_seq_quality_scores_data[2:]

    per_seq_quality_scores_df = pd.DataFrame([x.split("\t") for x in per_seq_quality_scores_data], columns=["Quality", "Count"], dtype=str)
    per_seq_quality_scores_df.to_csv(per_seq_quality_scores_data_path, sep="\t", index=False)



generate_general_stats_data(BEFORE_FASTQC, BEFORE_GENERAL_STATS)
generate_general_stats_data(AFTER_FASTQC, AFTER_GENERAL_STATS)

generate_per_seq_quality_scores_data(BEFORE_FASTQC, BEFORE_PER_SEQ_QUALITY_SCORES)
generate_per_seq_quality_scores_data(AFTER_FASTQC, AFTER_PER_SEQ_QUALITY_SCORES)