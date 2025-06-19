import os

# inputs 
INPUT_FASTQ = snakemake.input.raw_fastq

# outputs
OUT_FILE = snakemake.output.trimmed_fastq

# params
MINLEN_16S = snakemake.params.minlen_16S
MINLEN_18S = snakemake.params.minlen_18S
MINLEN_ITS = snakemake.params.minlen_ITS
HEADCROP = snakemake.params.headcrop
CROP = snakemake.params.crop
ORGANISM = snakemake.params.organism
GIVEN_THREAD = snakemake.params.t_thread


k = os.stat(INPUT_FASTQ)
k_size = k.st_size/(1024*1024)
print(k_size)
if k_size > 1000:
    avgqual =  "21"
elif k_size > 800:
    avgqual =  "20"
elif k_size > 700:
    avgqual =  "19"
elif k_size > 600:
    avgqual =  "18"
elif k_size > 500:
    avgqual =  "17"
elif k_size > 400:
    avgqual =  "16"
elif k_size > 300:
    avgqual =  "15"
elif k_size > 200:
    avgqual =  "13"
else:
    avgqual =  "11"


if ORGANISM == "Bacteria":
    command = f"trimmomatic SE -phred33 -threads {GIVEN_THREAD} {INPUT_FASTQ} {OUT_FILE} AVGQUAL:15 MINLEN:{MINLEN_16S} HEADCROP:{HEADCROP} CROP:{CROP}"
    os.system(command)

elif ORGANISM == "Fungi":
    command = f"trimmomatic SE -phred33 -threads {GIVEN_THREAD} {INPUT_FASTQ} {OUT_FILE} AVGQUAL:30 MINLEN:{MINLEN_ITS} HEADCROP:{HEADCROP} CROP:{CROP}"
    os.system(command)
else: 
    command = f"trimmomatic SE -phred33 -threads {GIVEN_THREAD} {INPUT_FASTQ} {OUT_FILE} AVGQUAL:{avgqual} MINLEN:{MINLEN_18S} HEADCROP:{HEADCROP} CROP:{CROP}"
    os.system(command)
