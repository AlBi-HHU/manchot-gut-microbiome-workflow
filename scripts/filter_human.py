import mappy as mp
from Bio.SeqIO.QualityIO import FastqGeneralIterator

a = mp.Aligner(snakemake.input['database'],preset='map-ont')  # load or build index
if not a: raise Exception("ERROR: failed to load/build index")

with open(snakemake.output['reads'],'w') as outfile,open(snakemake.output['stats'],'w') as statsfile:
    kept = 0
    filtered = 0
    for title, seq, qual in FastqGeneralIterator(snakemake.input['reads']): # read a fasta/q sequence
        try:
            hit = next(a.map(seq))
            filtered += 1
        except StopIteration:
            outfile.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))
            kept += 1
    statsfile.write('{}\t{}'.format(kept,filtered))

