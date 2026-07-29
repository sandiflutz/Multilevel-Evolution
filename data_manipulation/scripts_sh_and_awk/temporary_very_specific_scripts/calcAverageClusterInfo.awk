BEGIN{
	files=0;
	idt=0;
	maxtime=0.;
}{
	if(FNR==1){
		++files;
	}
# If not a comment, read record: 
	if ($1 !~ /#/){
		time[idt]=$1;
		if(time[idt]>maxtime)maxtime=time[idt];
		wcl[idt]+=$2;
		clsize[idt]+=$3;
		averw[idt]+=$4;
		nswcl[idt]+=$5;
		nclusters[idt]=$6;
		++sets[idt];
		++idt;
			
		if(idt>idtmax){
			idtmax=idt;
		}
	}
}END {

	printf("#1:time 2:<wb[ti]> 3:<averw[ti]> 4:<clsize[ti]> 5:<nswcl[ti]> 6:<nclusters[ti]> 7:sample\n");

	maxtime=time[idtmax-1];
	tinterval=int(idtmax/1000);
	if(tinterval==0)tinterval=1;
	timemeas0=maxtime-maxtime/10.0;
	tsample=0;
	meanwcl=0.0;
	meanwcl2=0.0;
	meanaverw=0.0;
	meanaverw2=0.0;
	meannswcl=0.0;
	meannclusters=0.0;
	savethis=0;
	for(i=0; i<idtmax; ++i){
		wcl[i]/=sets[i];
		averw[i]/=sets[i];
		clsize[i]/=sets[i];
		nswcl[i]/=sets[i];
		nclusters[i]/=sets[i];
	#	if(i==savethis){
		#	printf("%f %f %f %f %f %f %d\n",time[i],wcl[i],averw[i],clsize[i],nswcl[i],nclusters[i],sets[i]);
	#		savethis+=tinterval;
	#	}
		if(time[i]>timemeas0){
			meanwcl+=wcl[i];
			meanwcl2+=wcl[i]*wcl[i];
			meanaverw+=averw[i];
			meanaverw2+=averw[i]*averw[i];
			meannswcl+=nswcl[i];
			meannclusters+=nclusters[i];
			++tsample;
		}
	}

	stdwcl=(meanwcl2/tsample)-(meanwcl/tsample)*(meanwcl/tsample);		
	stdwcl=sqrt(stdwcl);
	stdaverw=(meanaverw2/tsample)-(meanaverw/tsample)*(meanaverw/tsample);		
	stdaverw=sqrt(stdaverw);
	#printf("\n #meanwcl=%f meanaverw=%f std(meanwcl)=%f std(meanaverw)=%f sample=%d\n",meanwcl/tsample,meanaverw/tsample,stdwcl,stdaverw,tsample);
	printf("%f %f %f %f %d\n",meanwcl/tsample,meanaverw/tsample,stdwcl,stdaverw,tsample);



}


