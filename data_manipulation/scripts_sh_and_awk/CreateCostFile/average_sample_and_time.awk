BEGIN{
	colunas=4;
	maximo=0;
}{	
	if(FNR==1){
	      ++files;
	}
# If not a comment, read record: 
	if(NF>colpars) colpars = NF;
	if(NF>colunas)colunas = NF;
	if($1 !~ /#/)	{
    		for(i=1; i <= NF ; ++i){
	 		coluna[imax]  += $i;
	 		coluna2[imax] += $i*$i;
         		++sets[imax]; 
	 		++imax;
     			if(imax>maximo)maximo = imax;
 		}
		if(files==1){
			for(i=11; i <=  NF; ++i){
				param[i-11]=$i;
			}
		}
	}

}END{
	

   #  printf("#1:nh/N 2:cost 3:mb 4:mh 5:rh 6:<w> 7:wb 8:clsizeb 9:std(<w>) 10:std(wb) 11:sample\n");
     imax=maximo/colunas;	
     
     timemax=coluna[(imax-1)*colunas]/sets[(imax-1)*colunas];
     timeM0=timemax-timemax/5.0;

     for(j=0; j<4; ++j){
	     aver[j]=0.;
	     aver2[j]=0.;
	     sample[j]=0;
     }

     for(i=0; i<imax; ++i){
	     for(j=0; j < 4; ++j){	
		     coluna[i*colunas+j]/=sets[i*colunas+j];
		     coluna2[i*colunas+j]/=sets[i*colunas+j];
		     if(coluna[i*colunas]>=timeM0){
		     	aver[j]+=coluna[i*colunas+j];
		     	++sample[j];
	     	     }
	     }
     }
     for(j=1; j < 4; ++j){	
	     aver[j]=aver[j]/(sample[j]);
     }
     printf("%f %f %f %f\n",param[1],aver[1],aver[2],aver[3]);
}


