from snakemake import shell

# inputs
FASTQ = snakemake.input.fastq

# outputs
FASTQC = snakemake.output.fastqc

mkdir_command = f"mkdir {FASTQC}"
fastqc_command = f"fastqc {FASTQ} --outdir {FASTQC}"
shell(mkdir_command)
shell(fastqc_command)
