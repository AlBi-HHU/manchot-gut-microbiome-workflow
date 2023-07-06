import glob
import pysam
import pandas as pd
import itertools

#https://alexwlchan.net/2018/12/iterating-in-fixed-size-chunks/
def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk

dfs = []
for identifiedElements,surroundingSequences in chunked_iterable(snakemake.input,2):

    samplename = identifiedElements.split('/')[-3]

    #print(patient,time)
    
    bestMappings={}
    associatedTaxa = {}
    
    
    with open(surroundingSequences,'r') as alignmentFile:
        next(alignmentFile)
        for line in alignmentFile:
            data = line.split('\t')
            identifier = data[0].split('_')
            sequence = identifier[0]
            adjacency = identifier[1]
            flags = "{0:012b}".format(int(data[1]))
            taxon = data[2]
            #9 unmapped
            #4 second in pair
            #3 not primary
            #2 read fails platform vendor checks
            if flags[2]=='1' or flags[3]=='1' or flags[4]=='1' or flags[9]=='1':
                continue
            else:
                if sequence not in bestMappings:
                    bestMappings[sequence] = {}
                if taxon not in bestMappings[sequence]:
                    bestMappings[sequence][taxon] = []
                bestMappings[sequence][taxon].append(adjacency)
                
        for sequence in bestMappings:
            associatedTaxa[sequence] = [
                taxon for taxon in bestMappings[sequence] if 
                all(adjacency in bestMappings[sequence][taxon] for adjacency in ['left','right'])
            ]
            
    tuples = []
    with open(identifiedElements,'r') as elementFile:
        for line in elementFile:
            data = line.split('\t')
            sequence = data[0]
            split = data[1].split('[')
            gene = split[0].strip()
            origin = split[1].strip()[:-1]
            other_annotations = split[2] if len(split) > 2 else 'None'
            assigned = False
            if sequence in associatedTaxa:
                for taxon in associatedTaxa[sequence]:
                    assigned = True
                    tuples.append((samplename,gene,origin,other_annotations,taxon,len(associatedTaxa[sequence])))
            else:
                    tuples.append((samplename,gene,origin,other_annotations,'Not Assigned',0))
    df = pd.DataFrame(tuples,columns=['Samplename','Gene','Source Taxon','Other Annotations','Assigned Taxon','Ambiguous Assignments'])
    dfs.append(df)
combined = pd.concat(dfs)
combined.to_csv(snakemake.output[0],index=False)
