awk '
BEGIN{
	split("A:T:C:G:N:del",bases,":"); 
	printf "chrom,pos,ref,major,minor";
}
(NR==1){ for(sampIX=4; sampIX<=NF; sampIX++) {printf ",af."$sampIX",rd."$sampIX} ; printf "\n" }
(NR>1){
	# record ref allele
	refBase=$3	
	
	# initialize baseTotals array for position
	for(ii=1;ii<=6;ii++){baseTotals[ii]=0}
	
	# store baseCts from each sample and add to totals
	for(ii=4;ii<=NF;ii++){
		split($ii,myArr,":");
		samp=ii-3
		for(base in myArr){
			baseCts[samp][base]=myArr[base]; 
			baseTotals[base]+=myArr[base];
		}
	} 
	
	# continue only if there are no deletions
	if(baseTotals[6]==0){
		
		# find major and minor alleles
		totaldepth=0.001; allelic=0
		minmaj[0]=""; minmaj[1]=""; 
		for(base=1;base<=4;base++) { totaldepth=totaldepth+baseTotals[base]} 
		for(base=1;base<=4;base++) { if((baseTotals[base]/totaldepth)>.01){ minmaj[allelic]=base; allelic++ } }
			
		# continue only if two alleles account for at least 99% of reads 
		if(allelic==2){
			# make major allele the one with more counts
			if(baseTotals[minmaj[0]]>baseTotals[minmaj[1]]){maj=minmaj[0]; min=minmaj[1]} else {maj=minmaj[1]; min=minmaj[0]}
			
			# continue only if the alt allele is at least 1% in two of the samples
			printMe=0; afrdString=""
			for(ii=4;ii<=NF;ii++){
				samp=ii-3
				sampRd=baseCts[samp][min]+baseCts[samp][maj]
				if(sampRd>0){ sampAF=baseCts[samp][min]/sampRd} else {sampAF=0}
				afrdString=afrdString","sampAF","sampRd				
				if(sampAF>=.01){printMe=printMe+1}
			}
						
			if(printMe>1){			
				# print minor allele freq for all samps	
				printf $1","$2","refBase","bases[maj]","bases[min] afrdString
				printf "\n"
			}
		}
	}
}'
			
		
	
