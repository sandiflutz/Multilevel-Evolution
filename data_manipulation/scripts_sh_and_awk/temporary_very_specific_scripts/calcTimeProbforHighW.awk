BEGIN{
	highw=0.85;
	files=0;
	bin=1.;
}{
	if(FNR==1){
		++files;
	}
# If not a comment, read record: 
	if ($1 !~ /#/){
		idt=int($1/bin);
		time[idt]=$1;
		wcl[idt]+=$2;
		clsize[idt]+=$3;
		averw[idt]+=$4;
		if(wcl[idt]>=highw){
			prob[idt]+=1.0;
		}
		++sets[idt];
			
		if(idt>idtmax){
			idtmax=idt;
		}
	}
}END {

	printf("#1:time 2:prob. of having a high investment cluster (wcl>%0.2f) 3:<wcl> 4:<averw> 5:<clsize> 6:sample\n",highw);

	for(i=0; i<idtmax; ++i){
		prob[i]/=sets[i];
		printf("%f %f %f %f %d %d\n",time[i],prob[i],wcl[i]/sets[i],averw[i]/sets[i],clsize[i]/sets[i],sets[i]);
	}



}


