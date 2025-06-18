import os


# inputs
INPUT_FASTQ = snakemake.input.input_fastq

# params
ORGANISM = snakemake.params.organism
OUTPUT_DIRECTORY = snakemake.params.outdir
GENOME_SIZE_BACTERIA = snakemake.params.genome_size_bacteria
GENOME_SIZE_VIRUS = snakemake.params.genome_size_virus
THREAD = snakemake.params.thread

if ORGANISM == "virus":
    genome_size = GENOME_SIZE_VIRUS
else:
    genome_size = GENOME_SIZE_BACTERIA

try:
    os.system(f"flye --nano-raw {INPUT_FASTQ} --threads 8 {genome_size} --asm-coverage 50 --out-dir {OUTPUT_DIRECTORY}")
except:
    os.system(f"flye --nano-raw {INPUT_FASTQ} --threads 8 --meta --out-dir {OUTPUT_DIRECTORY}")