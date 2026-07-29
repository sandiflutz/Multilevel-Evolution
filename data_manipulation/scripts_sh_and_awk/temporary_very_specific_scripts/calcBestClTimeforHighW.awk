BEGIN{
	idtmax=0;
	cp0=5;
	base[0]=0;
	idt=0;
	highw=0.85;

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

	printf("#1:initial time 2: final time 3:index 4:duration time of high cluster investment 5:mean average best cluster investment during the period 6:mean average system investment during the period 7:mean best cluster size during the period 8:time interval of the appearance of the first high investment cluster and the final time 9:maxdt\n");
	imax=idtmax-1;
	i=0;

	ndthw=0;
	dthighw[0]=0.0;
	meanaverw[0]=0.0;
	meanwcl[0]=0.0;
	averclsize[0]=0.0;
	nsample=0;
	maxdt=0.0;
	
	while(i<imax){
		if(wcl[i]>=highw){
			if((ndthw==0)&&(sample==0))time0[0]=time[i];
			timeF[ndthw]=time[i+1];

			dthighw[ndthw]+=time[i+1]-time[i];
			meanaverw[ndthw]+=averw[i];
			meanwcl[ndthw]+=wcl[i];
			averclsize[ndthw]+=clsize[i];
			++nsample;
			if(dthighw[ndthw]>maxdt)maxdt=dthighw[ndthw];
			if((wcl[i+1]<highw)||(i==imax-1)){
				if(i<imax-1)timeF[ndthw]=time[i];
				meanaverw[ndthw]/=nsample;
				meanwcl[ndthw]/=nsample;
				averclsize[ndthw]/=nsample;

				++ndthw;

				nsample=0;
				dthighw[ndthw]=0.0;
				meanaverw[ndthw]=0.0;
				meanwcl[ndthw]=0.0;
				averclsize[ndthw]=0.0;
				time0[ndthw]=time[i+1];
			}
		}
		++i;
	}

	for(i=0; i<ndthw; ++i){
		printf("%f %f %d %f %f %f %f %f %f\n",time0[i],timeF[i],i,dthighw[i],meanwcl[i],meanaverw[i],averclsize[i],time[imax]-time0[0],maxdt);
	}



}


