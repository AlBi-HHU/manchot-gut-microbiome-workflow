import pandas as pd

dfs = []

for f in snakemake.input:
	dfs.append(pd.read_csv(f))

pd.concat(dfs).to_csv(snakemake.output[0])