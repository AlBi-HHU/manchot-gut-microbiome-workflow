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
indexset = ["patientid","time"]
samples = pd.read_csv("samples.tsv",sep='\t', dtype={'patientid': str,'barcode':str,'time' : int}).set_index(indexset, drop=False)
validate(samples, "schemas/samples.schema.yaml")

# Verification
# Check for duplicates

dups = samples.index.duplicated( keep=False)
for s,d in zip(samples.index,dups):
    if d: #in this case we have a duplicate
        print("Duplicate entry for patientid: {} time: {} ... exiting!".format(s[0],s[1]))
        sys.exit(-1)
else:
    print("Sample sheet validated ... no duplicates!")

missingSamples = set()

outputList.append('data/output/sampleStats.csv')

#Kraken2/Bracken
if config['use_kraken_and_bracken']:
    outputList.append('data/output/mapping/KrakenFullDump.csv')
    outputList.append('data/output/mapping/KrakenFullDumpLangmead.csv')
    outputList+= expand('data/output/mapping/BrackenFullDump_{level}.csv',level=['P','G','S'])
#AMR    
if config['use_amr']:
    outputList.append('data/output/amr/fulldump.csv')

#Aggregated Subspecies Distances
if config['use_guttrsnp']:
    outputList.append('data/output/subspecies/distances.csv')
    outputList.append('data/output/subspecies/coverages.csv')

#Crassphage

if config['use_crassphage_verification']:
    outputList+=['data/output/crassphage/summary.csv']
    outputList+=['data/output/crassphage/mapping.csv']
    outputList+=['data/output/crassphage/migration.csv']

outputList+=['data/output/validation/combined.csv']

illumina_files = []



         
for idx,s in samples.iterrows():

    #Check if the required signals are available
    for _,sig in signals.iterrows():
        if sig['signalid'] == s['file'] or s['mode'] == 'fastq':
            if not sig['signalid'] in missingSignals:

                #Generate all output for samples that does not require additional illumina reads

                #Unassigned Reads Search
                if config['use_blast']:
                    outputList += ['data/auxiliary/unmappedReads/'+s['patientid']+'_'+str(s['time'])+'/blast.txt']
                    
                #Assembly
                if config['use_assembly']:
                    outputList+=['data/output/assembly/'+s['patientid']+'_'+str(s['time'])+'/prokka/annotation.gff']

                #MetaMaps
                if config['use_metamaps']:
                    outputList+=['data/output/kronaPlots/'+s['patientid']+'_'+str(s['time'])+'/metamaps.html']
                    
                #Zymo
                if s['patientid'].startswith('Zymo'):
                    outputList+=['data/auxiliary/mapping/zymo/'+s['patientid']+'_'+str(s['time'])+'.sam']                

                #outputList += ['data/auxiliary/samples/'+s['patientid']+'_'+str(s['time'])+'/reads.filtered.2.fastq']
                break
            else:
                missingSamples.add((s['patientid'],s['time']))
                print('We have to skip processing of patientid: {} time: {} because the signal file is missing (see above)'.format(s['patientid'],s['time']))
    else:
        print('Can\'t find the signal file {} required for sample patientid: {} time: {}, check the signal sheet!'.format(
        s['file'],s['patientid'],s['time']))
        sys.exit(-1)

    if s['illuminafile'] == 'none' or s['illuminafile'] == '' or s['illuminafile'] == 'None' or isinstance(s['illuminafile'],float):
        #This counts as no complementary illumina folder given in the sample sheet so there is no need for verification
        continue
    else:
        illumina_files.append(s)
        if config['use_markergene_alignments']:
            outputList+=['data/auxiliary/mapping/'+s['patientid']+'_'+str(s['time'])+'/markerGenes/'+s['patientid']+'_'+str(s['time'])+'.bam']
            if (
                os.path.isfile(os.path.join('data/input/shortreads',s['illuminafile']+'_R1.fastq.gz')) and
                os.path.isfile(os.path.join('data/input/shortreads',s['illuminafile']+'_R2.fastq.gz'))):
                continue
            else:
                print("Warning: Could not find the file {}_R[1 or 2].fastq.gz in 'data/input/shortreads' containing short reads for sample patientid: {} time: {}".format(s['illuminafile'],s['patientid'],s['time']))


for row in samples.itertuples():
    if row.mode == 'fastq':

        inputpath = 'data/input/reads/'+row.file +'/barcode'+str(row.barcode).zfill(2)


        #Sanity check
        if not os.path.exists(inputpath):
            print("Already demultiplexed read folder {} could not be found in data/input/reads but was specified in the samples.tsv table, aborting ...!".format(inputpath))
            sys.exit(-1)



### Rules

# Collect all results
rule all:
    input:
        outputList


### Include Additional Rules
include: "rules/common.snk"
include: "rules/basecalling.snk"
include: "rules/mapping.snk"
include: "rules/assembly.snk"
include: "rules/subspecies.snk"
include: "rules/verification.snk"
