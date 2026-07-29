BEGIN{
	idtmax=0;
	cp0=5;
	base[0]=0;
	}{
	if (FNR==1){
		idtmax=0;
		cp0=5;
		num=0;
		idt=0;
		base[0]=0;
	}
# If not a comment, read record: 
	if ($1 !~ /#/){
		time[idt]=$1;
		wcl[idt]=$2;
		clsize[idt]=$3;
		averw[idt]=$4;
		imax=clsize[idt];

		if(idt>0){
			base[idt]=base[idt-1]+clsize[idt-1];
		}
		
		id=base[idt]+num;
		idmax=base[idt]+clsize[idt];
		bestclist[id]=$(cp0+num);
		++num;
		if((base[idt]+num)>=idmax){
			num=0;
			++idt;
		}
		
		if(idt>idtmax){
			idtmax=idt;
		}
	}
}END {

	printf("#1:time 2:dt 3:corr[dt] 4:wcl[time] 5:wcl[time+dt] 6:averw[time] 7:averw[time+dt] 8:clsize[time] 9:clsize[time+dt]\n");
	corrt[0]=1.0;
	dtmax=200;
	imax=idtmax-dtmax;
	i=0;
	di=50;
	while(i<imax){
		for(dt=1; dt<dtmax; ++dt){
			j=i+dt;
			corr[dt]=0.0;
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
			printf("%f %f %f %f %f %f %f %d %d\n",i*5.0/100.0,dt*5.0/100.0,corrt[dt],wcl[i],wcl[i+dt],averw[i],averw[i+dt],clsize[i],clsize[i+dt]);
		}
		printf("\n");
		i+=di;
	}

}


