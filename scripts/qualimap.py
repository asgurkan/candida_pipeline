import os


bam=snakemake.input.bam

html=snakemake.params.html

k = os.stat(bam)
k_size = k.st_size/(1024*1024)

java_mem_size =  str(int(k_size/1024) + 1 )

java_mem_size1 = 32

java_mem_size_str = f" --java-mem-size={java_mem_size1}G"

command = f"flock -x data/qualimap.lock qualimap bamqc -bam {bam} --outdir {html} --outformat html {java_mem_size_str}"
#print(command)
os.system(command)
