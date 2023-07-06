from Bio.SeqIO.QualityIO import FastqGeneralIterator



with open(snakemake.output['reads'],'w') as outfile,open(snakemake.output['stats'],'w') as statsfile:
    kept = 0
    filtered = 0
    filter_read_ids = set()
    
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
            if taxid == '9606':
                filter_read_ids.add(read_id)
    
    
    for title, seq, qual in FastqGeneralIterator(snakemake.input['reads']): # read a fasta/q sequence
        title = title.split()[0]
        if title in filter_read_ids:
            filtered += 1
        else:
            outfile.write("@{}\n{}\n+\n{}\n".format(title, seq, qual))
            kept += 1
    statsfile.write('{}\t{}\t{}'.format(len(filter_read_ids),kept,filtered))

