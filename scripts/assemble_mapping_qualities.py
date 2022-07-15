import pysam
import pandas as pd

tuples = []

with open(snakemake.log[0],'w') as logfile:

    for f in snakemake.input:

        logfile.write('{}\n'.format(f))

        alignment_file = pysam.AlignmentFile(
            f
        )

        for aligned_segment in alignment_file:
            if not aligned_segment.is_unmapped:
                percentage_aligned = aligned_segment.query_alignment_length/ aligned_segment.query_length
                mq = aligned_segment.mapping_quality
                cigar_stats_nucleotide = aligned_segment.get_cigar_stats()[0]
                column_length = cigar_stats_nucleotide[0]+cigar_stats_nucleotide[1]+cigar_stats_nucleotide[2]
                mm_gaps = 1-(aligned_segment.get_tag('NM') / column_length)
                tuples.append(
                    (percentage_aligned,mq,mm_gaps)
                )

pd.DataFrame(tuples,columns=['Percentage Aligned','Mapping Quality','Identity']).to_csv(snakemake.output[0])
