import pandas as pd
import pysam
from scipy.stats import chi2


tuples = []
read_ids = set()
read_ids_unmapped = set()



al_file = snakemake.input['alignment']
kr_str_file = snakemake.input['stream']

with open(snakemake.log[0],'w') as logfile:


    unmapped = 0
    covered = {}
    filtered = 0
    identified = 0
    mapped = {}


    sample = al_file.split('/')[-1].rsplit('.',1)[0]

    alignment_file = pysam.AlignmentFile(
        al_file
    )


    for segment in alignment_file:

        if segment.is_unmapped:
            if segment.query_name not in read_ids_unmapped:
                read_ids_unmapped.add(segment.query_name)
                unmapped += 1
        else:

            #Check criteria for sufficiently sane alignment
            percentage_aligned = segment.query_alignment_length/ segment.query_length
            cigar_stats_nucleotide = segment.get_cigar_stats()[0]
            column_length = cigar_stats_nucleotide[0]+cigar_stats_nucleotide[1]+cigar_stats_nucleotide[2]
            mm_gaps = 1-(segment.get_tag('NM') / column_length) 

            if (mm_gaps < snakemake.config['mm2_perc_identity']) or (percentage_aligned < snakemake.config['mm2_perc_aligned']):
                if segment.query_name not in read_ids_unmapped:
                    read_ids_unmapped.add(segment.query_name)
                    filtered += 1
            else:
                
                if segment.query_name not in read_ids:
                    read_ids.add(segment.query_name)

                    if segment.mapping_quality < snakemake.config['mm2_uniqueness_mq']: #Low mapping quality, might not be unique to specific reference
                        identified += 1
                    else:

                        if not segment.reference_name in mapped:
                            mapped[segment.reference_name] = 0
                        mapped[segment.reference_name] += 1
    
                # Update Covered Positions anyways (Chimeric Reads)
                
                if segment.reference_name not in covered:
                    covered[segment.reference_name] = set()
                
                covered[segment.reference_name].update(
                    segment.get_reference_positions()
                )

    total = len(read_ids)+len(read_ids_unmapped)
    #del read_ids
    del read_ids_unmapped

    for ref in alignment_file.references:
        if ref not in covered:
            continue

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


        score = sum( (((dist-expected)**2) ) for dist in dists)/ expected

        chi_value = chi2.sf(score,1)

        covered_pos = len(covered[ref])

        ref_len = alignment_file.get_reference_length(ref)

        horizontal_coverage = covered_pos/ref_len
        
        if ref not in mapped: #only supplementary alignments that covered positions?
            continue
        
        mapped_reads = mapped[ref]
        
        fraction_reads = mapped_reads/total

        tuples.append(
            (
                sample,
                ref,
                horizontal_coverage,
                mapped_reads,
                fraction_reads,
                chi_value
            )
        )

    fraction_reads = identified/total
    tuples.append(
        (
            sample,
            'Ambiguous',
            None,
            identified,
            fraction_reads,
            None
        )        
    )   


    pd.DataFrame(tuples,columns=[
        'Sample','Reference','Horizontal Coverage','Mapped Reads','Fraction Mapped','Chi Values'
    ]).to_csv(snakemake.output['stats'])

    migration = {}

    with open(kr_str_file,'r') as infile:
        for l in infile:
            if l.startswith(('U','C')):
                data = l.split()
                readid = data[1]
                if readid in read_ids:
                    taxid = data[2]
                    if taxid not in migration:
                        migration[taxid] = 0
                    migration[taxid] += 1
                
    pd.DataFrame([(sample,k,v) for k,v in migration.items()],columns=['Sample','Tax ID','Reads']).to_csv(snakemake.output['migration'],index=False)
