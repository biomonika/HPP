#!/bin/bash
#SBATCH --job-name=extract_breakpoint_results.20240329
#SBATCH --partition=short
#SBATCH --mail-user=mcechova@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=1gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=extract_breakpoint_results.20240329.%j.log

#set -e
#set -x

pwd; hostname; date

mkdir -p T2T ass_stats mashmap

#copy T2T scaffold information
for a in cromwell-executions/findAssemblyBreakpoints/*/call-assessCompletness/execution/*.T2T.scaffolds.txt; do echo $a; b=`basename $a`; echo $b; cp -n $a T2T/$b; done;

#copy summary information
for a in cromwell-executions/findAssemblyBreakpoints/*/call-assessCompletness/execution/*.SUMMARY.txt; do echo $a; b=`basename $a`; echo $b; cp -n $a ass_stats/$b; done;

#copy mashmap files
for a in cromwell-executions/findAssemblyBreakpoints/*/call-evaluate/execution/*.mashmap.txt; do echo $a; b=`basename $a`; echo $b; cp -n $a mashmap/${b}; done;

#copy unified files
for a in cromwell-executions/findAssemblyBreakpoints/*/call-unifyAssembly/execution/*.unifiedAssembly.fa; do echo $a; b=`basename $a`; echo $b; cp -n $a unified/${b}; done;

#copy assembly statistics
for a in cromwell-executions/findAssemblyBreakpoints/*/call-assessCompletness/execution/*.T2T.scaffolds.txt; do b=`echo $a | sed "s/call-assessCompletness/call-createAssemblyStatistics/g" | sed "s/T2T.scaffolds/ass_stats/g"`; original=`echo $b | sed "s/\/[^/]*$/\/assembly.statistics.txt/g"`; new=`basename $b`; echo $original $new; cp -n $original ass_stats/${new}; done;

echo "Done."
date

