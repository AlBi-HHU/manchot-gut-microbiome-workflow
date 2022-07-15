import pandas as pd

pd.concat([pd.read_csv(x) for x in snakemake.input]).to_csv(snakemake.output[0])
