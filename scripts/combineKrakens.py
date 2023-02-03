import pandas as pd

tuples = []

for f in snakemake.input:
    
    with open(f,'r') as krakenfile:
        ll = krakenfile.read().splitlines()
        for l in ll:
            data = l.split('\t')
            readCount = int(data[1])
            level = data[3]
            taxonid = data[4]
            name = data[5].strip()
            fileinfo = f.split('/')
            iddata = fileinfo[-3].rsplit('_',1)
            patientid = iddata[0]
            time = int(iddata[1])
            tuples.append((patientid,time,readCount,level,name,taxonid))
    
df = pd.DataFrame(tuples,columns=['patientid','time','readcount','level','taxon','taxonid'])

df.to_csv(snakemake.output[0])

