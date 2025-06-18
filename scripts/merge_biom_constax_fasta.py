# inputs
INPUT_FILES = snakemake.input.fasta

# outputs
OUTPUT_FASTA = snakemake.output.fasta

main_read_dict = dict()

with open(INPUT_FILES) as f:
    data = f.readlines()
    for line in range(0,len(data),2):
        if data[line] not in main_read_dict:
            main_read_dict[data[line]] = data[line+1]
        else:
            pass
total = "".join([a+b for a,b in main_read_dict.items()])

with open(str(OUTPUT_FASTA),"w") as out:
    out.write(total)
