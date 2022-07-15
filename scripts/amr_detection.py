import pysam
from Bio import SeqIO

MIN_ALIGNMENT_LENGTH = snakemake.config['amr_min_alignment_length']
MIN_LENGTH_LR = snakemake.config['amr_min_length_sus']
MAX_NS_ALIGNED_SEGMENT = snakemake.config['amr_max_n_alignment']

#read database
reference_db = SeqIO.to_dict(
    SeqIO.parse(
        snakemake.input['card'],'fasta'
    )
)

f = pysam.AlignmentFile(snakemake.input['alignment'],'r')
secondary_supplementary_alignments = 0
short_alignments = 0
no_sufficient_surroundings = 0
unmapped_reads = 0
alignments_with_sus = 0
too_many_ns = 0
with open(snakemake.output['ide'],'w') as outfile, open(snakemake.output['sus'],'w') as susfile, open(snakemake.output['stats'],'w') as statsfile:
    segmentcounter = 0
    for aligned_segment in f:
        if aligned_segment.is_unmapped:
            unmapped_reads += 1
            continue
        if aligned_segment.is_secondary or aligned_segment.is_supplementary:
            secondary_supplementary_alignments += 1
            continue
        if aligned_segment.query_alignment_length < MIN_ALIGNMENT_LENGTH:
            short_alignments += 1
            continue
        if aligned_segment.query_alignment_sequence.count('N') >= MAX_NS_ALIGNED_SEGMENT:
            too_many_ns += 1
            continue
        else:
            left_sequence = aligned_segment.query_sequence[:aligned_segment.reference_start]
            right_sequence = aligned_segment.query_sequence[aligned_segment.reference_end:]
            lsl,rsl = len(left_sequence),len(right_sequence)
            if lsl < MIN_ALIGNMENT_LENGTH or rsl < MIN_ALIGNMENT_LENGTH:
                no_sufficient_surroundings += 1
            else:
                susfile.write('@{}_left\n'.format(segmentcounter))
                susfile.write('{}\n'.format(left_sequence))
                susfile.write('+\n')
                susfile.write('A'*len(left_sequence)+'\n')
                susfile.write('@{}_right\n'.format(segmentcounter))
                susfile.write('{}\n'.format(right_sequence))
                susfile.write('+\n')
                susfile.write('A'*len(right_sequence)+'\n')
                alignments_with_sus += 1
            #print(aligned_segment.reference_name)
            #resolve reference to full description
            full_reference = reference_db[aligned_segment.reference_name].description
            outfile.write('{}\t{}\n'.format(segmentcounter,full_reference))
            segmentcounter += 1
    statsfile.write('Unmapped Alignments\t{}\n'.format(unmapped_reads))
    statsfile.write('Secondary or Supplementary Alignments\t{}\n'.format(secondary_supplementary_alignments))
    statsfile.write('Short Alignments\t{}\n'.format(short_alignments))
    statsfile.write('Alignments Without Sufficient Surrounding Sequence Length\t{}\n'.format(no_sufficient_surroundings))
    statsfile.write('Alignments With Sufficient Surrounding Sequence Length\t{}\n'.format(alignments_with_sus))
    statsfile.write('Alignments With Excessive Amount of Ns\t{}\n'.format(too_many_ns))
