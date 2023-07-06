from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import random
from functools import lru_cache
import pandas as pd
import os
import shutil

with open(snakemake.log[0],'w') as logfile:

    #clean incomplete runs
    if os.path.exists(snakemake.output['reads']):
        logfile.write('removing incomplete run ...')
        shutil.rmtree(snakemake.output['reads'])
        logfile.write('done!')
        logfile.flush()
    os.makedirs(snakemake.output['reads'],exist_ok=True)

    #Load Taxonomy

    idtonames = {}

    names = snakemake.input['taxonomy_names']

    nodes = snakemake.input['taxonomy_nodes']


    with open(names,'r') as f:
        for l in f.read().splitlines():
            d=[x.strip() for x in l.split('|')]
            if d[3] == 'scientific name':
                idtonames[d[0]] = d[1]


    taxonomy = {}


    levels = {}

    with open(nodes,'r') as f:
        for l in f.read().splitlines():
            d=[x.strip() for x in l.split('|')]
            taxonomy[d[0]] = d[1]
            levels[d[0]] = d[2]



    idtonames['0'] = 'Unassigned'

    RANK_NAMES = [ "root",
                   "superkingdom",
                   "kingdom",
                   "subkingdom",
                   "superphylum",
                   "phylum",
                   "subphylum",
                   "superclass",
                   "class",
                   "subclass",
                   "infraclass",
                   "superorder",
                   "order",
                   "suborder",
                   "infraorder",
                   "parvorder",
                   "superfamily",
                   "family",
                   "subfamily",
                   "tribe",
                   "subtribe",
                   "genus",
                   "subgenus",
                   "species group",
                   "species subgroup",
                   "species",
                   "subspecies",
                   "varietas",
                   "forma",
                   "strain"]

    @lru_cache(maxsize=128)
    def resolve_highest_ancestor_below_or_at(current_node,last_valid_node=None,level='genus'):
        current_level = levels[current_node]
        
        #ranks that are not in our list are not evaluated we treat them as "no rank"
        current_level = 'no rank' if current_level not in RANK_NAMES else current_level
        
        # if we observe a valid rank below our target level we memorize it as a potential valid node to fall back
        if current_level != 'no rank' and RANK_NAMES.index(current_level) > RANK_NAMES.index(level):
            last_valid_node = current_node
            
            
        # if the level is the target level we return it
        if current_level == level:
            return current_node
        # if the level is above the target level we need the last node below genus
        elif current_level != 'no rank' and RANK_NAMES.index(current_level) < RANK_NAMES.index(level):
            return last_valid_node
        elif current_node == '1':
            return last_valid_node
        else:
            return resolve_highest_ancestor_below_or_at(
                taxonomy[current_node],
                last_valid_node,
                level
            )

    logfile.write('Generating Verification Map ...\n')
    logfile.flush()

    verification_map = {}

    with open(snakemake.input['stream'],'r') as kraken_stream:
        next(kraken_stream)
        for line in kraken_stream:
            if not line.startswith(('U\t','C\t')):
                continue
            #print(line)
            data = line.split('\t')
            read_id = data[1]
            #Discard unassigned reads
            if data[0] == 'U':
                continue
            #Resolve corresponding genus or highest assigned node below genus
            taxid = data[2]
            ancestor = resolve_highest_ancestor_below_or_at(taxid,level=snakemake.wildcards['level'])
            verification_map[read_id] = ancestor

    logfile.write('Done ... Generating suitable read set \n')
    logfile.flush()
    
    random.seed(160922)

    reads =  FastqGeneralIterator(open(snakemake.input['reads']))
    suitable_reads = set()

    for title, seq, qual in reads:
        sid = title.split()[0]
        if sid not in verification_map:
            continue
        else:
            suitable_reads.add((title,seq,qual))
            
    logfile.write('Done, sampling ...\n')
    logfile.flush()

    limit = min(100000,len(suitable_reads))

    read_selection = random.sample(sorted(suitable_reads),limit)
    
    logfile.write('Done, writing reads \n')
    logfile.flush()
        

    for title, seq, qual in read_selection:
        sid = title.split()[0]

        ancestor_id = verification_map[sid]

        with open(snakemake.output['reads']+'/{}.fq'.format(ancestor_id),'a') as outfile:
            outfile.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))
            
