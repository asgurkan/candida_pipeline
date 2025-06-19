input = snakemake.input
output = snakemake.output

new_fasta_content = []
with open(str(input)) as f:
    for line in f:
        for s_line in line.split("\t"):
            new_fasta_content.append(s_line.rstrip())
with open(str(output), 'w') as filehandle:
    filehandle.writelines("%s\n" % place for place in new_fasta_content)
