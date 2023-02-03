import pandas as pd
import pysam
from scipy.stats import chi2


tuples = []
ref_lens = {}


with open(snakemake.log[0],'w') as logfile:

    
    for f in snakemake.input:
        
        unmapped = 0
        covered = {}
        filtered = 0

        
        taxon = f.split('/')[-1].split('.')[0].split('_')[-1]

        logfile.write('{}\n'.format(f))

        sample = f.split('/')[-1].rsplit('.',1)[0]

        alignment_file = pysam.AlignmentFile(
            f
        )
        
        processed_read_ids = set()
        
        for segment in alignment_file:
                
            if segment.is_unmapped:
                if segment.query_name not in processed_read_ids:
                    processed_read_ids.add(segment.query_name)
                    unmapped += 1
            else:
                
                #Check criteria for sufficiently sane alignment
                percentage_aligned = segment.query_alignment_length/ segment.query_length
                cigar_stats_nucleotide = segment.get_cigar_stats()[0]
                column_length = cigar_stats_nucleotide[0]+cigar_stats_nucleotide[1]+cigar_stats_nucleotide[2]
                mm_gaps = 1-(segment.get_tag('NM') / column_length) 
                
                if (mm_gaps < snakemake.config['mm2_perc_identity']) or (percentage_aligned < snakemake.config['mm2_perc_aligned']):
                    if segment.query_name not in processed_read_ids:
                        processed_read_ids.add(segment.query_name)
                        filtered += 1

                else:
                    if segment.reference_name not in covered:
                        covered[segment.reference_name] = set()
                
                    covered[segment.reference_name].update(
                        segment.get_reference_positions()
                    )
                    
                if segment.query_name not in processed_read_ids:
                    processed_read_ids.add(segment.query_name)


        total = len(processed_read_ids)
                
        covered_pos = sum(len(x) for x in covered.values())
        
        ref_len = sum(alignment_file.get_reference_length(ref) for ref in alignment_file.references)
                
        horizontal_coverage = covered_pos/ref_len
        mapped_reads = 1-(
            filtered+unmapped
        )/total
            
        chi_values = []
            
        for ref in alignment_file.references:
            if ref not in covered:
                continue

            # chi squared test ; 0 Hypothese: Abweichung ist Random
            # Score bedeutet: Wahrscheinlichkeit solche oder eine noch stärkere Abweichung zu sehen, 
            # unter Annahme der Gleichverteilung

            cur = 0
            dist = 1
            dists = []

            while cur < alignment_file.get_reference_length(ref):
                if cur in covered[ref]:
                    dists.append(dist)
                    dist = 1
                else:
                    dist += 1
                cur += 1

            dists.append(dist)


            number_of_hits = len(dists)
            expected = sum(dists)/number_of_hits
            #print(expected)

            score = sum( (((dist-expected)**2) ) for dist in dists)/ expected
            #print(score)
            chi_values.append(chi2.sf(score,1))

        tuples.append(
            (
                sample,
                taxon,
                horizontal_coverage,
                mapped_reads,
                chi_values,
                filtered,
                unmapped,
                total
            )
        )
        


pd.DataFrame(tuples,columns=[
    'Sample','Taxon','Horizontal Coverage','Mapped Reads','Chi Values','Filtered Reads','Unmapped Reads','Total Reads'
]).to_csv(snakemake.output[0],index=False)

