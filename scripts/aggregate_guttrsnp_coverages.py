import pandas as pd

tbls = []

for f in snakemake.input:
    df = pd.read_csv(f,sep='\t')
    split1 = f.split('/')
    idstring = split1[-1].replace('.pileup.summary','')

    df['Samplename'] = idstring
    tbls.append(df)
pd.concat(tbls).to_csv(snakemake.output[0],index=False)
