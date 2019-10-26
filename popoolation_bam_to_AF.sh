#!/bin/bash

ml biology; ml samtools; ml java; 
bamDirs=$1
afDir=$2
ref=${3:-dmel}
refFasta=$($(dirname $0)/../get_ref_fasta.sh $ref)
chrom=$4
threads=$5
encoding=${6:-sanger}
scriptDir=${7:-$(dirname $0)}

popoolationDir=$PI_HOME/SOFTWARE/popoolation2_1201


bamfiles=$(echo $bamDirs | tr ',' '\n' | sed 's/$/\/*.dedup.realign*bam/')
echo "calculating raw allele frequencies for:"
echo $bamfiles | tr ' ' '\n'

headerText=$(echo "chrom pos ref "$(for bam in $bamfiles; do basename $bam | sed 's/.dedup.*.bam//'; done)  | tr ' ' '\t')

java -ea -Xmx7g -jar $popoolationDir/mpileup2sync.jar \
--fastq-type $encoding --threads $threads \
--input <(samtools mpileup -f $refFasta -r $chrom --min-MQ 10 --min-BQ 20 -R $bamfiles ) \
--output >( (echo $headerText && cat) | $scriptDir/parse_sync_to_AF.sh > $afDir/allSamps.$chrom.af) 

### get average read depth for each samp 
ncol=$(head -1 $afDir/allSamps.$chrom.af | tr ',' ' ' | wc -w)
paste 	<(echo $bamfiles | tr ' ' '\n')  \
	<(tail -n +2 $afDir/allSamps.$chrom.af | cut -f$(echo `seq 7 2 "$ncol"` | tr ' ' ',') -d',' | awk -F ',' '
		{for(ii=1;ii<=NF;ii++){cts[ii]=cts[ii]+$ii}} 
		END {for(samps in cts){print samps"\t"cts[samps]/NR}}
	' | sort -k1g) \
> $afDir/allSamps.$chrom.af.meanrd

