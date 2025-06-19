import pandas as pd

shannon = snakemake.input.shannon 
chao1=snakemake.input.chao1
simpson=snakemake.input.simpson
chao1_ci=snakemake.input.chao1_ci
ENSPIE=snakemake.input.ENSPIE
fisher=snakemake.input.fisher
mcintosh_d =snakemake.input.mcintosh_d
mcintosh_e =snakemake.input.mcintosh_e
output = snakemake.output


shannon_df = pd.read_csv(str(shannon), sep="\t",index_col=0)
chao1_df = pd.read_csv(str(chao1), sep="\t",index_col=0)
simpson_df = pd.read_csv(str(simpson), sep="\t",index_col=0)
chao1_ci_df=pd.read_csv(str(chao1_ci), sep="\t",index_col=0)
ENSPIE_df=pd.read_csv(str(ENSPIE), sep="\t",index_col=0)
fisher_df =pd.read_csv(str(fisher), sep="\t",index_col=0)
mcintosh_d_df=pd.read_csv(str(mcintosh_d), sep="\t",index_col=0)
mcintosh_e_df=pd.read_csv(str(mcintosh_e), sep="\t",index_col=0)
result = pd.concat([shannon_df, chao1_df, chao1_ci_df, simpson_df, ENSPIE_df, fisher_df, mcintosh_d_df, mcintosh_e_df], axis=1)
result.to_csv(str(output), sep="\t")


    

