from Bio import SeqIO
import pandas as pd
import numpy as np

readLengths = []

for record in SeqIO.parse(snakemake.input[0], 'fastq'):
	readLengths.append(len(record))

pd.DataFrame(
	[(snakemake.wildcards.samplename,
	  len(readLengths),
	  np.sum(readLengths),
	  np.median(readLengths),
	  np.mean(readLengths),
	  np.std(readLengths),
	  np.min(readLengths),
	  np.max(readLengths)
	)],
	columns=['samplename', 'nReads', 'nBases','median','mean','standard deviation','minimum','maximum']
).to_csv(snakemake.output[0])
