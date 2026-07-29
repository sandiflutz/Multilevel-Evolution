BEGIN{
	idtmax=0;
	cp0=5;
	base[0]=0;
	idt=0;
	}
	function ceil(x,   y) { y=int(x); return(x>y ? y+1 : y) }{
# If not a comment, read record: 
	if ($1 !~ /#/){
		for(i=1; i<=NF; ++i){
			if(i==1){
				time[idt]=$i;
			}else if(i==2){
				wcl[idt]=$i;
			}else if(i==3){
				clsize[idt]=$i;
				if(idt>0){
					base[idt]=base[idt-1]+clsize[idt-1];
				}
			}else if(i==4){
				averw[idt]=$i;
			}else{
				id=base[idt]+i-cp0;
				bestclist[id]=$i;
			}
		}
		++idt;
		
		if(idt>idtmax){
			idtmax=idt;
		}
	}
}END {

	printf("#1:time 2:frac. hosts of best cluster that remain from the last time step 3:wcl 4:averw 5:clsize[last time] 6:clsize[current time]\n");
	i=0;
	while(i<idtmax-1){
		j=i+1;
		frac[j]=0.;
		for(k1=0; k1<clsize[i]; ++k1){
			for(k2=0; k2<clsize[j]; ++k2){
				if(bestclist[base[i]+k1]==bestclist[base[j]+k2]){
					frac[j]+=1.0;
				}
			}
		}
		frac[j]/=clsize[i];
		
		printf("%f %f %f %f %d %d\n",j*5.0/100.0,frac[j],wcl[j],averw[j],clsize[i],clsize[j]);
		++i;
	}

}


