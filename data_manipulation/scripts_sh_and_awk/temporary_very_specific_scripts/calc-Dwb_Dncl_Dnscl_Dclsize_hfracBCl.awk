BEGIN{
	idt=0;
}
function fabs(x){
	if(x<0.){
		return (-1.0)*x;
	}else{
		return x;
	}
}{
# If not a comment, read record: 
        if ($1 !~ /#/){
                time[idt]=$1;
                wcl[idt]=$2;
                clsize[idt]=$3;
                averw[idt]=$4;
		nsimilarw[idt]=$5;
		ncluster[idt]=$6;
		if(idt>0){
			base[idt]=base[idt-1]+clsize[idt-1];
                }
		for(i=7; i<=NF; ++i){
			id=base[idt]+i-7;
			bestclist[id]=$i;
		}


                if(idt>idtmax){
                        idtmax=idt;
                }
		++idt;
        }
}END {

        printf("#1:time 2:wb[i]-wb[i-1] 3:<w>[i]-<w>[i-1] 4:clsize[i]-clsize[i-1] 5:nsimilaw[i]-nsimilarw[i-1] 6:ncluster[i]-ncluster[i-1] 7:hosts frac. that remain in the best cluster 8:wb 9:<w> 10:maxdiffnsw\n");

	frac[0]=0.0;
	maxdiffnsw=1;
	for(i=0; i<idtmax-1; ++i){
		j=i+1;
		frac[j]=0.;
		for(ki=0; ki<clsize[i]; ++ki){
			for(kj=0; kj<clsize[j]; ++kj){
				if(bestclist[base[i]+ki]==bestclist[base[j]+kj]){
					frac[j]+=1.0;
				}
			}
			frac[j]/=clsize[i];
		}
		if((fabs(nsimilarw[i+1]-nsimilarw[i])>maxdiffnsw)&&(time[i+1]>=1000.))maxdiffnsw=fabs(nsimilarw[i+1]-nsimilarw[i]);
	}

        for(i=1; i<idtmax; ++i){
                printf("%f %f %f %d %d %d %f %f %f %f\n",time[i],wcl[i]-wcl[i-1],averw[i]-averw[i-1],clsize[i]-clsize[i-1],nsimilarw[i]-nsimilarw[i-1],ncluster[i]-ncluster[i-1],frac[i],wcl[i],averw[i],maxdiffnsw);
        }
}

