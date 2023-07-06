import pysam
import glob
import itertools as it
from Bio import SeqIO
import os
import pandas as pd

hc_mapping = {}

for f in glob.glob(snakemake.input['reads']+'/*.fq'): #Iterate over 
    
    ancestor_id = f.split('/')[-1].split('.')[0] 
    
    reads = SeqIO.parse(f,'fastq')

    if ancestor_id == 'None':
        for read in reads:
            if 'uv' not in hc_mapping:
                hc_mapping['uv'] = 0
            hc_mapping['uv'] += 1  
    else:
        if not os.path.exists(snakemake.params['alignments']+'/{}.bam'.format(ancestor_id)):
            for read in reads:
                if ancestor_id not in hc_mapping:
                    hc_mapping[ancestor_id] = 0
                hc_mapping[ancestor_id] += 1
        else:
            
            data = {}
            
            alignment_file = pysam.AlignmentFile(
                snakemake.params['alignments']+'/{}.bam'.format(ancestor_id),check_sq=False
            )
            
            for segment in alignment_file:
                if segment.query_length < 500:
                    continue

                if segment.query_name not in data:
                    data[segment.query_name] = {
                        'intervals' : [],
                        'identities' : [],
                        'read_length' : 0,
                        'unmapped' : False
                    }
                    
                if segment.is_unmapped:
                    data[segment.query_name]['unmapped'] = True
                    continue
                    
                data[segment.query_name]['read_length'] = segment.query_length

                q_st = segment.query_alignment_start
                q_en = segment.query_alignment_end
                data[segment.query_name]['intervals'].append((q_st,q_en))
                cigar_stats_nucleotide = segment.get_cigar_stats()[0]
                column_length = cigar_stats_nucleotide[0]+cigar_stats_nucleotide[1]+cigar_stats_nucleotide[2]
                mm_gaps = 1-(segment.get_tag('NM') / column_length) 
                data[segment.query_name]['identities'].append(mm_gaps)

                    
            for read in data:
                if data[read]['unmapped']:
                    if ancestor_id not in hc_mapping:
                        hc_mapping[ancestor_id] = (0,0)
                    hc_mapping[ancestor_id] = (hc_mapping[ancestor_id][0],hc_mapping[ancestor_id][1]+1)                    
                else:
                    if (len(data[read]['intervals']) > 1):
                        #check for overlaps
                        change_detected = True
                        while(change_detected):
                            change_detected = False
                            for iv1,iv2 in it.combinations(range(len(data[read]['intervals'])),2):
                                iv1_st,iv1_en = data[read]['intervals'][iv1]
                                iv2_st,iv2_en = data[read]['intervals'][iv2]
                                if ((iv1_st > iv2_st and iv1_en < iv2_en) or #within
                                   (iv1_st < iv2_st and iv1_en > iv2_en) or #around
                                    (iv1_st < iv2_st and iv1_en > iv2_st) or #left overlap
                                    (iv1_st < iv2_en and iv1_en > iv2_en) #right overlap
                                   ):
                                    #print('Found overlap {} / {}'.format(intervals[iv1],intervals[iv2]))
                                    #expand
                                    data[read]['intervals'][iv1] = (
                                        min(iv1_st,iv2_st),
                                        max(iv1_en,iv2_en)
                                    )        
                                    del data[read]['intervals'][iv2]
                                    change_detected=True
                                    break
                    covered = 0
                    for interval in data[read]['intervals']:
                        covered += (interval[1]-interval[0])+1
                    covered /= data[read]['read_length']
                    average_identity = sum(data[read]['identities'])/len(data[read]['identities'])

                    if ancestor_id not in hc_mapping:
                        hc_mapping[ancestor_id] = (0,0)
                    if covered >= snakemake.config['mm2_perc_aligned'] and average_identity >= snakemake.config['mm2_perc_identity']:
                        hc_mapping[ancestor_id] = (hc_mapping[ancestor_id][0]+1,hc_mapping[ancestor_id][1])
                    else:
                        hc_mapping[ancestor_id] = (hc_mapping[ancestor_id][0],hc_mapping[ancestor_id][1]+1)

                        
result_tuples = []

for k,v in hc_mapping.items():
    if k == 'uv':
        continue
    result_tuples.append((k,v))

results_refs = pd.DataFrame(result_tuples,columns=['Taxon ID','Validation'])
results_refs

def calculate_rate(row):
    if isinstance(row['Validation'],int):
        return
    total = 0
    validated = 0

    validated += row['Validation'][0]
    total += row['Validation'][0]+row['Validation'][1]
    return validated/total

applied = results_refs.apply(calculate_rate,axis=1)
print(applied)
results_refs['Validation Rate'] = applied
results_refs = results_refs.sort_values(by=['Validation Rate'])
results_refs.to_csv(snakemake.output['aggregation'],index=False)
