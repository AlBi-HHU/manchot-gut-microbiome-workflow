#!/usr/bin/env bash

while getopts "r:o:f:k:c:g:e:" opt
do
  case $opt in
    r) rootfolder="$OPTARG" ;;
    o) outfolder="$OPTARG" ;;
    f) flowcell="$OPTARG" ;;
    k) kit="$OPTARG" ;;
    c) cpus="$OPTARG" ;;
    g) gpu="$OPTARG" ;;
    e) exec="$OPTARG" ;;
  esac
done


folder=$(find -L "$rootfolder"/ -name 'final_summary*.txt' -printf '%h\n') && \
echo "Identified $folder as the root folder to start basecalling" && \

resumeArgument=""

if ls "${outfolder}"/*.log 1> /dev/null 2>&1; then
    echo "A previous basecalling run was found, attempting to resume basecalling"
    resumeArgument="--resume"
else
    echo "This will be the first attempt at basecalling for this sample!"
fi && \

if [ "$gpu" -eq 1 ]
then
  ${exec} --recur --input_path ${folder}/fast5_pass/ --save_path $outfolder --flowcell $flowcell --kit $kit --device auto $resumeArgument
else
  ${exec} --recur --input_path ${folder}/fast5_pass/ --save_path $outfolder --flowcell $flowcell --kit $kit --cpu_threads_per_caller $cpus $resumeArgument
fi
