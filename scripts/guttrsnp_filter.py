import pysam
import sys
import cython

@cython.cfunc
def main(
    alignment_file,
    output_file,
    vertical_cutoff,
    horizontal_cutoff,
    threads,
    log_file
    ):
    
    total_refs : cython.int
    ref_len : cython.int
    ref_tid : cython.int
    with open(output_file,'w') as outfile, open(log_file,'w') as logfile:
    
        logfile.write('beginning decompression ...')
        logfile.flush()
        f = pysam.AlignmentFile(alignment_file,'rb',threads=threads,check_sq=False)

        total_refs = len(f.references)
        logfile.write('Total of {} references detected\n'.format(total_refs))
        logfile.flush()
        for contig in f.references:

            ret = []
        
            ref_tid = f.get_tid(contig)
            
            #logfile.write(str(ref_tid)+'/'+str(total_refs))
            ref_name = f.get_reference_name(ref_tid)
            
            #potential speedup
            reads = f.count(ref_name)
            if reads == 0:
                continue
                
            ref_len = f.get_reference_length(contig)

            cnt_a,cnt_c,cnt_g,cnt_t = f.count_coverage(ref_name)
            if len(cnt_a) != ref_len:
                logfile.write('Coverage Array does not match the length of contig')
                logfile.flush()
                sys.exit(-1)
                
            pos : cython.int
            
            for pos,alleles in enumerate(zip(cnt_a,cnt_c,cnt_g,cnt_t)):
                a,c,g,t = alleles
                if a+c+g+t < vertical_cutoff:
                    continue
                else:
                    ret.append('{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            ref_name,
                            pos,
                            a,
                            c,
                            g,
                            t
                        )
                    )
                    
            if len(ret)/ref_len < horizontal_cutoff:
                logfile.write('skipping reference {} due to insufficient horizontal coverage ... \n'.format(ref_name))
                logfile.flush()
                continue
            else:
                for entry in ret:
                    outfile.write(entry)

        logfile.write('Script finished!')   
    
    
    
if "snakemake" in locals():
    vertical_cutoff = snakemake.config['guttrsnp_filter_vertical']
    horizontal_cutoff = snakemake.config['guttrsnp_filter_horizontal']
    main(
        snakemake.input["alignment"],
        snakemake.output["alleles"],
        vertical_cutoff,
        horizontal_cutoff,
        int(snakemake.resources['cpus']),
        snakemake.log[0]
    )
else:
    main(
        sys.argv[1],
        sys.argv[2],
        int(sys.argv[3]),
        float(sys.argv[4]),
        int(sys.argv[5]),
        sys.argv[6]
    )
