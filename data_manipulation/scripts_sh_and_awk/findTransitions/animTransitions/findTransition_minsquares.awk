BEGIN{
	imax=0;
}
function fabs(x){ return(x<0.0 ? -x : x)}
function ceil(x,   y) { y=int(x); return(x>y ? y+1 : y) }
{
	
	if($1 i!~ /#/){
		for(i=1; i<NF; ++i){
			if(i==1){
				x[imax]=$i;
			}else if(i==2){
				y[imax]=$2;
				++imax;
			}
		}
	}
	if(imax > lines)lines=imax;

}END{
	if(PRINT_AxT!=0){
		printf("#1:t0 2:tmax 3:a 4:b 5:ymed 6:stdymed 7:errymed 8:dadt 9:samplewindow\n");
	}
	if(samplewindow==0){
		expo=ceil((log(x[lines-1])/log(10)));
		mul=ceil((x[lines-1]/10^(expo-2))-10);
		samplewindow=int((10+mul)*10^(expo-2)/100);
		PRINT_TTRANS=1;
		PRINT_SLOPE=0;
	}else{
		PRINT_TTRANS=0;
		PRINT_SLOPE=1;
	}
	sampleinterval=1;
	
	imax=lines-samplewindow;
	j0=0;
	jmax=samplewindow;
	bins=0;
	while(j0<imax){
		xmed=0.0;
		ymed=0.0;
		ymed2=0.0;
		for(j=j0; j<jmax; ++j){
			xmed+=x[j];
			ymed+=y[j];
			ymed2+=y[j]*y[j];
		}
		xmed/=samplewindow;
		ymed/=samplewindow;
		ymed2/=samplewindow;
		fact=0.0;
		facb=0.0;
		for(j=j0; j<jmax; ++j){
			fact+=x[j]*(y[j]-ymed);
			facb+=x[j]*(x[j]-xmed);
		}
		a[bins]=fact/facb;
		b[bins]=ymed-a[bins]*xmed;
		errymed=0.0;
		for(j=j0; j<jmax; ++j){
			errymed+=(y[j]-(a[bins]*x[j]+b[bins]))*(y[j]-(a[bins]*x[j]+b[bins]));
		}
		errymed/=samplewindow;
		errymed=sqrt(errymed);

		stdymed=ymed2-ymed*ymed;
		stdymed=sqrt(stdymed);
		if(PRINT_AxT!=0){
			printf("%f %f %0.16f %f %f %f %0.16f %d\n",x[j0],x[jmax-1],a[bins],b[bins],ymed,stdymed,errymed,samplewindow);
		}
		j0+=sampleinterval;
		jmax=j0+samplewindow;
		++bins;
	}
	
	i=0;
	ok=0;
	amax=0.;
	amin=1.;
	for(i=0; i<bins; ++i){
		id0=i+i*(sampleinterval-1);
		idf=i+i*(sampleinterval-1)+samplewindow;
		if(a[i]>amax){
			amax=a[i];
			b_amax=b[i];
			x0_amax=x[id0];
			y0_amax=y[id0];
		}
		if(i>0){
			aabs=fabs(a[i]);
			if((aabs<amin)&&(a[i-1]>0.0)&&(a[i+1]<0.0)&&(ok==0)){
				amin=aabs;
				b_amin=b[i];
				x0_amin=x[id0];
				y0_amin=y[id0];
				id0_amin=id0;
				ok=1;
			}
		}

	}
	nsample=0;
	yasyn=0.0;
	yasyn2=0.0;
	for(i=id0_amin; i<lines; ++i){
		yasyn+=y[i];
		yasyn2+=y[i]*y[i];
		++nsample;
	}
	yasyn/=nsample;
	yasyn2/=nsample;
	stdy=yasyn2-yasyn*yasyn;
	stdy=sqrt(stdy);

	if(PRINT_AMAX!=0){
		if(PRINT_SLOPE!=0){
			printf("%f %f %d\n",amax,b_amax,samplewindow,file_index);
		}else{
			printf("%f %f %f %f %f %f %f %f %d %d\n",x0_amax,x0_amin,amax,b_amax,y0_amax,y0_amin,yasyn,stdy,samplewindow,file_index);
		}
	}
}
