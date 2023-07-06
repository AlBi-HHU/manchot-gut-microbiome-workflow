import pandas as pd

tbls = []

for f in snakemake.input:
    df = pd.read_csv(f)
    split1 = f.split('/')
    idstring = split1[-1].replace('.dist','')
    df['Distance Name'] = idstring
    tbls.append(df)
pd.concat(tbls).to_csv(snakemake.output[0],index=False)
