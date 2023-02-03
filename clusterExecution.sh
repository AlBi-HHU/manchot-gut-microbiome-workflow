#!/bin/bash

umask 0002

module load Miniconda/3

conda activate snakemake

unset PYTHONHOME
unset PYTHONPATH

while getopts "p:j:" opt
do
  case $opt in
    p) projectID="$OPTARG" ;;
    j) maxNrOfConcurrentJobs="$OPTARG" ;;
  esac
done

if [ -z "$projectID" ] || [ -z "$maxNrOfConcurrentJobs" ];
then
  echo "Usage: clusterExecution -p projectID -j maxNrOfConcurrentJobs"
  exit 1
fi

mkdir -p clusterLogs

type snakemake >/dev/null 2>&1 || { echo >&2 "I require snakemake but it's not installed or added to your path.  Aborting..."; exit 1; }

snakemake --jobs $maxNrOfConcurrentJobs --rerun-incomplete --printshellcmds --use-conda --cluster-status "python cluster/statuscommand.py" --reason --jobscript cluster/jobscript.sh --cluster "qsub -e clusterLogs/{rule}.{wildcards}.{jobid}.errors -o clusterLogs/{rule}.{wildcards}.{jobid}.output -A ${projectID} -l select=1:ncpus={resources.cpus}:ngpus={resources.gpus}:mem={resources.mem} -l walltime={resources.walltime}"

