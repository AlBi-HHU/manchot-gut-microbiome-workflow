import pandas as pd

tbls = []

for f in snakemake.input:
    df = pd.read_csv(f)
    split1 = f.split('/')
    idstring = split1[-1].replace('.dist','')
    split3 = idstring.split('_')
    sourceid,sourcetime = split3[0],split3[1]
    destid,desttime = split3[3],split3[4]    
    df[['Source ID','Source Time','Destination ID','Destination Time']] = sourceid,sourcetime,destid,desttime
    tbls.append(df)
pd.concat(tbls).to_csv(snakemake.output[0],index=False)
