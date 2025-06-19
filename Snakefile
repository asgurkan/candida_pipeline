import os
from glob import glob
import pandas as pd
import numpy as np
import itertools

configfile: "config.yaml"

c_auris_ref_fna = "/home/massbiome/Desktop/pipelines/candida_pipeline/workspace/ref_genomes1/GCF_003013715.1_ASM301371v2_genomic.fna"
c_auris_ref_gff = "/home/massbiome/Desktop/pipelines/candida_pipeline/workspace/ref_genomes1/genbank.genbank"
dbs_folder = '/home/massbiome/Desktop/pipelines/candida_pipeline/dbs'
blast_custom_clade_db = "/home/massbiome/Desktop/pipelines/candida_pipeline/dbs/blast_c_auris_clade_db/c_auris_clades_db"

samples = pd.read_csv("samples.tsv", sep="\t")

def get_samples():
    """Returns list of samples."""
    return list(samples['sample'].unique())
    
def get_fastq(wildcards):
    sample_name = str(wildcards.sample)
    fastqs = samples.query("sample == @sample_name")[["fastq_path"]].iloc[0]
    ret_str = f"{fastqs.fastq_path}"
    return ret_str

def get_projects_with_samples(w_project_sample):
    proje_sample = []
    project_name = w_project_sample
    samples_df = samples.query(f"proje=='{project_name}'")

    for a_key,a_value in samples_df.iterrows():
        proje_sample.append(a_value["sample"])

    print(proje_sample)
    return proje_sample

rule all:
    input:
        #expand(["data/quality_check/before_trim_quality_graph_data/{sample}_per_seq_quality_scores_before_trim.tsv"], sample=get_samples())
        #expand(["data/wgs/denovo_alignment/{sample}/assembly.fasta"], sample=get_samples())
        #expand(["data/gambit/{sample}_results.csv"], sample=get_samples())
        expand(["data/wgs/blast/filtered/{sample}_blast_result.txt"], sample=get_samples())
        #expand(["data/samples/{sample}.fastq"], sample=get_samples())

rule merge_first_fastqs:
   input:
       get_fastq
   conda:
        "envs/general.yml"
   output:
       fastq="data/samples/{sample}.fastq"
   shell:
       """
       cat {input}/*.fastq > {output.fastq}
       """

rule fastq_quality_check_before_trim:
    input:
        fastq = "data/samples/{sample}.fastq"
    output:
        fastqc = directory("data/quality_check/{sample}_before_trim")
    params:
        sample_name = lambda wildcards: wildcards.sample
    conda:
       "envs/general.yml"
    script:
        "scripts/check_fastqc.py"
    

rule trimmomatic:
   input:
       raw_fastq = "data/samples/{sample}.fastq"
   output:
       trimmed_fastq = "data/trimmed/{sample}.fastq"
   params:
       minlen_16S =config["trimmomatic"]["minlen_16S"], 
       minlen_18S =config["trimmomatic"]["minlen_18S"],
       minlen_ITS =config["trimmomatic"]["minlen_ITS"],
       headcrop =config["trimmomatic"]["headcrop"], 
       crop = config["trimmomatic"]["crop"], 
       t_thread = config["thread"],
       organism = config["organism"]
   conda:
       "envs/general.yml"
   message: 
        "Trimmomatic has been started!"
   script:
       "scripts/trimmomatic.py"


rule fastq_quality_check_after_trim:
    input:
        fastq = "data/trimmed/{sample}.fastq"
    output:
        fastqc = directory("data/quality_check/{sample}_after_trim")
    params:
        sample_name = lambda wildcards: wildcards.sample
    conda:
       "envs/general.yml"
    script:
        "scripts/check_fastqc.py"


rule create_quality_graph_data:
    input:
        before_fastqc = directory("data/quality_check/{sample}_before_trim"),
        after_fastqc = directory("data/quality_check/{sample}_after_trim")
    output:
        before_per_seq_quality_scores = "data/quality_check/before_trim_quality_graph_data/{sample}_per_seq_quality_scores_before_trim.tsv",
        before_general_stats = "data/quality_check/before_trim_quality_graph_data/{sample}_general_stats.tsv",
        after_per_seq_quality_scores = "data/quality_check/after_trim_quality_graph_data/{sample}_per_seq_quality_scores_before_trim.tsv",
        after_general_stats = "data/quality_check/after_trim_quality_graph_data/{sample}_general_stats.tsv"
    params:
        sample_name = lambda wildcards: wildcards.sample
    conda:
       "envs/general.yml"
    script:
        "scripts/create_quality_graph_data.py"

rule de_novo_assembly: 
    input: 
        input_fastq = "data/samples/{sample}.fastq"
    params:
        organism = config["organism"],
        thread = config["thread"],
        outdir = "data/wgs/denovo_alignment/{sample}",
        genome_size = config["candida_genome_size"]
    output: 
        de_novo_fasta = "data/wgs/denovo_alignment/{sample}/assembly.fasta"
    conda: 
        srcdir("envs/flye.yml")    
    shell:
        "flye --nano-raw {input.input_fastq} --threads {params.thread} {params.genome_size} --out-dir {params.outdir}"

rule blast_custom_db:
    input: 
        assembly_file = "data/wgs/denovo_alignment/{sample}/assembly.fasta"
    output:
        blast_result = "data/wgs/blast/raw/{sample}_blast_result.txt"
    params:
        blast_db = f"{dbs_folder}/blast_c_auris_clade_db/c_auris_clades_db"
    conda:
        srcdir("envs/blast.yml")   
    shell:
        """
        blastn -query {input.assembly_file} -db {params.blast_db} -outfmt 6 -out {output.blast_result} -culling_limit 1
        """

rule blast_processing:
    input:
        blast_result = "data/wgs/blast/raw/{sample}_blast_result.txt"
    output:
        combined_result = "data/wgs/blast/filtered/{sample}_blast_result.txt"
    params:
        dbs_folder = dbs_folder
    conda:
        srcdir("envs/general.yml")
    script:
        "scripts/blast_processing.py"
        

# rule gambit:
#     input:
#         assembly_fasta = "data/wgs/denovo_alignment/{sample}/assembly.fasta"
#     output:
#         gambit_result_csv = "data/gambit/gambit_results.csv"
#     params:
#         gambit_db = config["gambit_db"]
#     conda:
#        "envs/gambit.yml"
#     shell:
#         "gambit -d {params.gambit_db} query -o {output.gambit_result_csv} {input.assembly_fasta}"

# rule minimap2:
#     input:
#         fastq = "data/samples/{sample}.fastq",
#         reference = c_auris_ref_fna
#     params:
#         preset = "-x map-ont", # Oxford Nanopore to reference mapping (-k15)
#         threads = config["thread"],
#     output:
#         sam = "data/wgs/alignment/{sample}/{sample}.sam",
#         bam = "data/wgs/alignment/{sample}/{sample}.bam",
#         sorted_bam = "data/wgs/alignment/{sample}/{sample}.sorted.bam",
#         sorted_bai = "data/wgs/alignment/{sample}/{sample}.sorted.bam.bai",
#     conda: 
#         srcdir("envs/bcftools.yml")
#     shell:
#         """
#         # Index the reference genome
#         samtools faidx {input.reference}

#         # Run minimap2 to align reads and generate SAM file
#         minimap2 -a --frag=yes -t {params.threads} {params.preset} {input.reference} {input.fastq} > {output.sam}

#         # Convert SAM to BAM (without read group)
#         samtools view -bS {output.sam} > {output.bam}

#         # Add read group using samtools addreplacerg
#         samtools addreplacerg -r "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}" -o {output.bam}.rg {output.bam}

#         # Sort the BAM file
#         samtools sort {output.bam}.rg -o {output.sorted_bam}

#         # Index the sorted BAM file
#         samtools index {output.sorted_bam}
#         ""


def get_best_clade_reference(sample):
    blast_file = f"data/wgs/blast/{sample}_blast_result.txt"
    blast_results = pd.read_csv(blast_file, sep=' ', header=None)
    blast_results.columns = ['Clade', '% Identity', 'Path']
    blast_results['% Identity'] = blast_results['% Identity'].astype(float)
    best_clade_row = blast_results.loc[blast_results['% Identity'].idxmax()]
    return best_clade_row['Path']



rule minimap2:
    input:
        fastq = "data/samples/{sample}.fastq",
        reference = lambda wildcards: get_best_clade_reference(wildcards.sample)
    params:
        preset = "-x map-ont", 
        threads = config["thread"],
    output:
        sam = "data/wgs/alignment/{sample}/{sample}.sam",
        bam = "data/wgs/alignment/{sample}/{sample}.bam",
        sorted_bam = "data/wgs/alignment/{sample}/{sample}.sorted.bam",
        sorted_bai = "data/wgs/alignment/{sample}/{sample}.sorted.bam.bai",
    conda: 
        srcdir("envs/bcftools.yml")
    shell:
        """
        samtools faidx {input.reference}
        minimap2 -a --frag=yes -t {params.threads} {params.preset} {input.reference} {input.fastq} > {output.sam}
        samtools view -bS {output.sam} > {output.bam}

        # Add read group using samtools addreplacerg
        samtools addreplacerg -r "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}" -o {output.bam}.rg {output.bam}

        samtools sort {output.bam}.rg -o {output.sorted_bam}
        samtools index {output.sorted_bam}
        """

rule qualimap:
    input:
        bam = "data/wgs/alignment/{sample}/{sample}.sorted.bam"
    output:
        directory("data/wgs/alignment/qualimap_results/{sample}_qualimap_result")
    conda: 
        srcdir("envs/qualimap.yml")
    params:
        html="data/wgs/alignment/qualimap_results/{sample}_qualimap_result"
    script:
       "scripts/qualimap.py"

rule variant_calling:
    input:
        # qualimap = rules.qualimap.output,
        reference_fna = c_auris_ref_fna, 
        reference_gff = c_auris_ref_gff, 
        index = "data/wgs/alignment/{sample}/{sample}.sorted.bam.bai",
        align = "data/wgs/alignment/{sample}/{sample}.sorted.bam" 
    output: 
        snippy_out_dir = directory("data/wgs/variants/{sample}/")
    params:
        prefix = "{sample}_snps"
    threads: 
        config["thread"]
    conda:
        srcdir("envs/snippy.yml")
    shell:
        "snippy --outdir {output.snippy_out_dir} --ref {input.reference_gff} --bam {input.align}"
