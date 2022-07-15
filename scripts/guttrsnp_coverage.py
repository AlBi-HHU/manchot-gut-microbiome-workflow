import pandas as pd

metadata = pd.read_csv(snakemake.input["marker_genes_metadata"],dtype={'Taxon ID' : str}).set_index('Taxon ID')

output = {}

with open(snakemake.input['alleles'],'r') as infile:

    for l in infile:
        reference,pos,a,c,g,t = l.split()
        taxon,gene = reference.split('.',maxsplit=1)
        if not taxon in output:
            output[taxon] = []
        output[taxon].append(
            sum(int(x) for x in [a,c,g,t])
        )

with open(snakemake.output['summary'],'w') as outfile: 

    outfile.write(
        '{}\t{}\t{}\n'.format(
            'Taxon ID','Average Vertical Coverage','Horizontal Coverage'
        )  
    )
    for taxon in output:
        total_length = metadata.loc[taxon]['Total Length']
        outfile.write('{}\t{}\t{}\n'.format(
                taxon,sum(output[taxon])/len(output[taxon]),len(output[taxon])/total_length
            )
        )
