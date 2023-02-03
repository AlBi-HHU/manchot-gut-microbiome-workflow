import sys
import pandas as pd

def main(
    metadata_marker_genes,
    pileup_source_path,
    pileup_dest_path,
    output_path,
    vertical_cutoff_abs,
    vertical_cutoff_rel_source,
    vertical_cutoff_rel_dest
    ):
    
    metadata = pd.read_csv(metadata_marker_genes,dtype={'Taxon ID' : str}).set_index('Taxon ID')
    
    print('parsing source ...')
    #parse source
    src = {}
    with open(pileup_source_path,'r') as infile:
        for l in infile:
            gene,pos,a,c,g,t = l.split()
            a = int(a)
            c = int(c)
            g = int(g)
            t = int(t) 
            total = a+c+g+t
            #track present alleles
            src[(gene,pos)] = (
                a/total >= vertical_cutoff_rel_source,
                c/total >= vertical_cutoff_rel_source,
                g/total >= vertical_cutoff_rel_source,
                t/total >= vertical_cutoff_rel_source,
            )
    print('done')
    results = {}

    with open(pileup_dest_path,'r') as infile:
        for l in infile:
            gene,pos,a,c,g,t = l.split()
            taxon = gene.split('.')[0]
            
            if taxon not in results:
                results[taxon] = {
                    'overlapping_present_positions' : 0,
                    'high_confidence_alleles' : 0,
                    'missing_src_positions' : 0,
                    'available_src_positions' : 0,
                    'likely_not_present' : 0,
                    'likely_present' : 0
                }
            
            a = int(a)
            c = int(c)
            g = int(g)
            t = int(t)      
            results[taxon]['overlapping_present_positions'] += (gene,pos) in src
            total = a+c+g+t
            alleles = [a,c,g,t]
            for idx,allele in enumerate(alleles):
                #determine high confidence allele set
                if (allele > vertical_cutoff_abs) and (allele/total > vertical_cutoff_rel_dest):
                    results[taxon]['high_confidence_alleles'] += 1
                    if (gene,pos) in src:
                        results[taxon]['available_src_positions'] += 1
                        if src[(gene,pos)][idx]:
                            results[taxon]['likely_present'] += 1
                        else:
                            results[taxon]['likely_not_present'] += 1
                    else:
                        results[taxon]['missing_src_positions'] += 1
                        
        outframe = pd.DataFrame(
            ((
                taxon,
                results[taxon]['high_confidence_alleles'],
                results[taxon]['missing_src_positions'],
                results[taxon]['available_src_positions'],
                results[taxon]['likely_not_present'],
                results[taxon]['likely_present'],
                results[taxon]['overlapping_present_positions'],
                results[taxon]['overlapping_present_positions']/metadata.loc[taxon]['Total Length'],
                results[taxon]['likely_not_present']/results[taxon]['overlapping_present_positions'] if results[taxon]['overlapping_present_positions'] != 0 else 'NaN'
            ) for taxon in results) , columns=[
                'Taxon',
                'High Confidence Alleles',
                'Missing Source Positions',
                'Available Source Positions',
                'Likely Not Present Alleles',
                'Likely Present Alleles',
                'Overlapping Present Positions',
                'Share Of Overlap',
                'GutTrSnp Distance'
            ]
 
        )                  
        outframe.to_csv(output_path,index=False)
                                                                                  
                                                                                             
if "snakemake" in locals():
    vertical_cutoff_abs = snakemake.params['VerticalCutoffAbsolute']
    vertical_cutoff_rel_source = snakemake.params['VerticalCutoffRelativeSource']
    vertical_cutoff_rel_dest = snakemake.params['VerticalCutoffRelativeDestination']
    main(
        snakemake.input["marker_genes_metadata"],
        snakemake.input["pileup_source"],
        snakemake.input["pileup_dest"],
        snakemake.output["distance"],
        vertical_cutoff_abs,
        vertical_cutoff_rel_source,
        vertical_cutoff_rel_dest
    )
else:
    main(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3],
        sys.argv[4],
        int(sys.argv[5]),
        float(sys.argv[6]),
        float(sys.argv[7])
    )
