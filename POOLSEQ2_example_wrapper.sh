#!/bin/bash

ml system rclone biology samtools
rootDir=$1
ref=$2
gdrive=$3
scriptDir=${4:-$GROUP_HOME/poolseq_pipelines}
if [ ! -e $rootDir/rawAF ]; then mkdir -p $rootDir/rawAF; fi
chroms=(2L 2R 3L 3R X)

for chrom in ${chroms[*]}; do

	time="20:00:00"; mem=10000; threads=12
	cmd="$scriptDir/PIPELINE2a_bam_to_rawAF/popoolation_bam_to_AF.sh $rootDir/mapped $rootDir/rawAF $ref $chrom $threads; rclone copy $rootDir/rawAF $gdrive/rawAF --include=allSamps.$chrom.af"; 
	$scriptDir/submit_sbatch_jobs.sh "$cmd" $chrom rawAF $time $mem 1 $threads none $rootDir/jobScripts
done
