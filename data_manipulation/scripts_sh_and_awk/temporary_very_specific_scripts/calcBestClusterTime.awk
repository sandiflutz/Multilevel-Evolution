BEGIN{
	idtmax=0;
	cp0=5;
	base[0]=0;
	idt=0;
	}
	function ceil(x,   y) { y=int(x); return(x>y ? y+1 : y) }{
# If not a comment, read record: 
	if ($1 !~ /#/){
		time[idt]=$1;
                wcl[idt]=$2;
                clsize[idt]=$3;
                averw[idt]=$4;
                nswcl[idt]=$5;
                nclusters[idt]=$6;
                if(idt>0){
                        base[idt]=base[idt-1]+clsize[idt-1];
                }
                for(i=7; i<=NF; ++i){
                        id=base[idt]+i-cp0;
                        bestclist[id]=$i;
                }
		++idt;
		
		if(idt>idtmax){
			idtmax=idt;
		}
	}
}END {

	printf("#1:time 2:duration time of best cluster 3:wcl[initial time] 4:wcl[final time] 5:averw[initial time] 6:averw[final time] 7:clsize[initial time] 8:clsize[final time] 9:mean time correlation\n");
	corrt[0]=1.0;
	imax=idtmax-10;
	i=0;

	while(i<imax){
		dt=0;
		meancorr=0.;
		do{
			++dt;
			j=i+dt;
			#time correlation calculation
			if(clsize[i]>clsize[j]){
				maxclsize=clsize[i];
			}else{
				maxclsize=clsize[j];
			}
			corrt[dt]=0.0;
			for(k1=0; k1<clsize[i]; ++k1){
				for(k2=0; k2<clsize[j]; ++k2){
					if(bestclist[base[i]+k1]==bestclist[base[j]+k2]){
						corrt[dt]+=1.0;
					}
				}
			}
			corrt[dt]/=maxclsize;
			if(corrt[dt]>0.)meancorr+=corrt[dt];
		}while((corrt[dt]>0.)&&(i+dt<idtmax));
		
		if(meancorr>0.)meancorr/=(dt-1);
		if(corr[dt]==0.){
			printf("%f %f %f %f %f %f %d %d %f\n",i*5.0/100.0,(dt-1.0)*5.0/100.0,wcl[i],wcl[i+dt-1],averw[i],averw[i+dt-1],clsize[i],clsize[i+dt-1],meancorr);
		}
		i+=dt;
	}

}


