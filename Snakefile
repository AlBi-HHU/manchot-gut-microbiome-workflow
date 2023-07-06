import pandas as pd
import os
import sys
from glob import glob

from snakemake.utils import validate

### Load Configuration Files

print('Loading and validating configuration files ...')

configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")

### This list will contain all files that are desired
outputList = []


print('Configuration files validated successfully!')

### Load Signal Sheet

signals = pd.read_csv("signals.tsv",sep='\t').set_index("signalid", drop=False)
validate(signals, "schemas/signals.schema.yaml")

missingSignals = []
#Check if all signals given exist as either folders or symlinks
for idx,s in signals.iterrows():
    if os.path.exists(os.path.join('data/input/signals',s['signalid'])) or os.path.exists(os.path.join('data/input/signals',s['signalid']+'.tar.gz')):
        pass
    else:
        missingSignals.append(s['signalid'])
        print("Warning: Could not find the folder/file {} or {}.tar.gz in 'data/input/signals' containing fast5 data.".format(s['signalid'],s['signalid']))


# Read Sample Sheet

samples = pd.read_csv("samples.tsv",sep='\t', dtype={'samplename': str,'barcode':str}).set_index('samplename', drop=False)
validate(samples, "schemas/samples.schema.yaml")

# Verification
# Check for duplicates

dups = samples.index.duplicated( keep=False)
for s,d in zip(samples.index,dups):
    if d: #in this case we have a duplicate
        print("Duplicate entry for sample: {} ... exiting!".format(s))
        sys.exit(-1)
else:
    print("Sample sheet validated ... no duplicates!")

missingSamples = set()



#Kraken2/Bracken
if config['use_kraken_and_bracken']:
    outputList.append('data/output/mapping/'+config['project_prefix']+'_KrakenFullDump.csv')
    outputList+= expand('data/output/validation/'+config['project_prefix']+'_{level}.csv',level=['species','genus'])
    #outputList+= expand('data/output/mapping/'+config['project_prefix']+'_BrackenFullDump_{level}.csv',level=['P','G','S'])
#AMR    
if config['use_amr']:
    outputList.append('data/output/amr/'+config['project_prefix']+'_fulldump_amr.csv')

#VF
if config['use_vf']:
    outputList.append('data/output/vf/'+config['project_prefix']+'_fulldump_vf.csv')

#Aggregated Subspecies Distances
if config['use_guttrsnp']:
    outputList.append('data/output/subspecies/'+config['project_prefix']+'_distances.csv')
    outputList.append('data/output/subspecies/'+config['project_prefix']+'_coverages.csv')

#Crassphage

if config['use_crassphage_verification']:
    outputList+=['data/output/crassphage/summary.csv']
    outputList+=['data/output/crassphage/mapping.csv']
    outputList+=['data/output/crassphage/migration.csv']

#Sample Statistics Pre/Post Filtering

outputList.append('data/output/'+config['project_prefix']+'_sampleStats.reads.csv')
outputList.append('data/output/'+config['project_prefix']+'_sampleStats.reads.filtered.2.csv')

illumina_files = []


for idx,s in samples.iterrows():
    add_files = False
    if s['mode'] == 'fastq':
        add_files = True #ToDo: Prüfe ob, fastq existiert
    else:     
    #Check if the required signals are available
        for _,sig in signals.iterrows():
            if sig['signalid'] == s['file']:
                if (sig['signalid'] not in missingSignals):
                    add_files = True
                    break
                else:
                    missingSamples.add(s['samplename'])
                    print('We have to skip processing of samplename: {} because the signal file is missing (see above)'.format(s['samplename']))
                    sys.exit(-1)
        else:
            print('Can\'t find the signal file {} required for samplename: {}, check the signal sheet!'.format(
            s['file'],s['samplename']))
            sys.exit(-1)

    if add_files == True:
        #Generate all output for samples that does not require additional illumina reads

        #Unassigned Reads Search
        if config['use_blast']:
            outputList += ['data/auxiliary/unmappedReads/'+s['samplename']+'/blast.txt']

        #Assembly
        if config['use_assembly']:
            outputList+=['data/output/assembly/'+s['samplename']+'/prokka/annotation.gff']

        #MetaMaps
        if config['use_metamaps']:
            outputList+=['data/output/kronaPlots/'+s['samplename']+'/metamaps.html']

        #Zymo
        if s['samplename'].startswith('Zymo'):
            outputList+=['data/auxiliary/mapping/zymo/'+s['samplename']+'.sam']   
            
        #Gene-Set Coverage Lists
        if config['use_vf']:
            outputList+=['data/auxiliary/mapping/'+s['samplename']+'/vf/coverage.txt']           
        if config['use_amr']:
            outputList+=['data/auxiliary/mapping/'+s['samplename']+'/amr/coverage.txt'] 
            
        #Illumina-Based Analysis
        if s['illuminafile'] == 'none' or s['illuminafile'] == '' or s['illuminafile'] == 'None' or isinstance(s['illuminafile'],float):
            #This counts as no complementary illumina folder given in the sample sheet so there is no need for verification
            continue
        else:
            illumina_files.append(s)
            if config['use_markergene_alignments']:
                outputList+=['data/auxiliary/mapping/'+s['samplename']+'/markerGenes/'+s['samplename']+'.bam']
                if (
                    os.path.isfile(os.path.join('data/input/shortreads',s['illuminafile']+'_R1.fastq.gz')) and
                    os.path.isfile(os.path.join('data/input/shortreads',s['illuminafile']+'_R2.fastq.gz'))):
                    continue
                else:
                    print("Warning: Could not find the file {}_R[1 or 2].fastq.gz in 'data/input/shortreads' containing short reads for samplename: {}".format(s['illuminafile'],s['samplename']))

#Keep track of barcode notation mode, depeneding on sequencing platform the folders can be named barcodexx or xx
barcode_notation_mode = {}

for row in samples.itertuples():
    if row.mode == 'fastq':

        path_option_a = 'data/input/reads/'+row.file +'/barcode'+str(row.barcode).zfill(2)
        path_option_b = 'data/input/reads/'+row.file +'/'+str(row.barcode).zfill(2)

        #Sanity check
        if not (os.path.exists(path_option_a) or os.path.exists(path_option_b)):
            print("Neither {} nor {} could be found but were specified in the samples.tsv table!".format(path_option_a,path_option_b))
        else:
            if os.path.exists(path_option_a):
                barcode_notation_mode[row.samplename] = glob(path_option_a+'/*.fastq')+glob(path_option_a+'/*.fastq.gz')
            else:
                barcode_notation_mode[row.samplename] = glob(path_option_b+'/*.fastq')+glob(path_option_b+'/*.fastq.gz')           

### Rules

# Collect all results
rule all:
    input:
        outputList


### Include Additional Rules
include: os.path.join('rules','common.snk')
include: os.path.join("rules","basecalling.snk")
include: os.path.join("rules","mapping.snk")
include: os.path.join("rules","assembly.snk")
include: os.path.join("rules","subspecies.snk")
include: os.path.join("rules","verification.snk")
