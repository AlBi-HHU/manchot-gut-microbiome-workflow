import glob
from multiprocessing import Pool
import os

os.makedirs(snakemake.params['alignments'],exist_ok=True)
os.makedirs(snakemake.params['prefix'],exist_ok=True)
def align(f):
    
    ancestor_id = f.split('/')[-1].split('.')[0]
    
    if ancestor_id == 'None':
        return

    if not os.path.exists(snakemake.input['database']+'/{}.mmi'.format(ancestor_id)):
        print('no database for id: {}'.format(ancestor_id))
        return
    
    if not os.path.exists(snakemake.params['alignments']+'/{}.bam'.format(ancestor_id)):
        addonargs = ''
        #if ancestor_id in ['9605','5506']:
        addonargs = ' --split-prefix '+snakemake.params['prefix']+'/'+ancestor_id
        print('mapping for id:{}'.format(ancestor_id))

        minimap_cmd = 'minimap2 -a -x map-ont -M 0 -t 1 --hard-mask-level {} --secondary=no {}/{}.mmi {} | samtools sort -o {}/{}.bam'.format(addonargs,snakemake.input['database'],ancestor_id,f,snakemake.params['alignments'],ancestor_id)

        return_value = os.system(minimap_cmd)
        
        if return_value != 0:
            raise Exception('Minimap CMD failed')
        
with Pool(int(snakemake.resources['cpus'])) as pool:
    pool.map(align,glob.glob(snakemake.input['reads']+'/*.fq'))
    
with open(snakemake.output['flag'],'w') as outfile:
    outfile.write('done!')
