from glob import glob
from Bio import SeqIO

import os

#Assemble valid eukaryota refs

eukaryota_ids = []

for f in glob(os.path.join(snakemake.input['database'],'*.mmi')):
    eukaryota_ids.append(f.split('/')[-1].split('.')[0])


#Phase 1: Extract 'S' level hits with > 100 Reads

kraken_report = snakemake.input['report']
kraken_stream = snakemake.input['stream']

tax_ids_to_check = set()

with open(kraken_report,'r') as infile:
    for l in infile:
        d = l.split('\t')
        if (d[3] == 'S') and (d[4] in eukaryota_ids):
            reads = int(d[2])
            if reads >= 100:
                tax_ids_to_check.add(d[4])
        else:
            continue

relevant_reads = {}

with open(kraken_stream,'r') as infile:
    for l in infile:
        d = l.split('\t')
        if d[0] == 'C':
            if d[2] in tax_ids_to_check:
                relevant_reads[d[1]]=d[2]
        else:
            continue
            
del tax_ids_to_check
      
reads = SeqIO.parse(
    snakemake.input['reads'],'fastq'
)

reads_sorted = {}


for read in reads:
    if read.id in relevant_reads:
        if relevant_reads[read.id] not in reads_sorted:
            reads_sorted[relevant_reads[read.id]] = []
        reads_sorted[relevant_reads[read.id]].append(read)

os.makedirs(
    snakemake.output['eukaryotes'],
    exist_ok=True
)

for taxon in reads_sorted:

    with open(
        os.path.join(
            snakemake.output['eukaryotes'],
            '{}.fastq'.format(taxon)
        ),'w'
    ) as outfile:
        SeqIO.write(
            reads_sorted[taxon],outfile,'fastq'
        )

