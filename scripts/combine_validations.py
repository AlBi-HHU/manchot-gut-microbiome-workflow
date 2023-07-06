import pandas as pd
import numpy as np

#Load Taxonomy

idtonames = {}

with open(snakemake.input['names'],'r') as f:
    for l in f.read().splitlines():
        d=[x.strip() for x in l.split('|')]
        if d[3] == 'scientific name':
            idtonames[d[0]] = d[1]

levels = {}
        
with open(snakemake.input['nodes'],'r') as f:
    for l in f.read().splitlines():
        d=[x.strip() for x in l.split('|')]
        levels[d[0]] = d[2]

tables = []

for f in snakemake.input['validations']:
    sample = f.split('/')[-2]
    data = pd.read_csv(f,dtype={'Taxon ID' : str,'Validation' :str})
    data['Sample'] = sample
    tables.append(data)
combined = pd.concat(tables)
def label(x):
    ret = x
    if x in idtonames:
        ret += ' '+idtonames[x]
    if levels[x] != 'genus':
        ret += '(Pseudogenus)'
    return ret
combined['Taxon Name'] = combined['Taxon ID'].apply(label)
combined['Reads'] = combined['Validation'].apply(
    lambda x : eval(x)[0]+eval(x)[1] if isinstance(eval(x),tuple) else int(x)
)
combined['Exponent'] = np.log10(combined['Reads']).round().astype(int)
totals = combined.groupby('Sample')['Reads'].sum()
combined['Fraction Estimate'] = combined.apply(
    lambda x : x['Reads'] / totals[x['Sample']],axis=1
)
combined.to_csv(snakemake.output['combined'],index=False)
