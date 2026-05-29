BEGIN{
	imax=0;
	}{
	if($1 !~ /#/){
		if(colunas<NF)colunas=NF;
		for(i=1; i<=NF; ++i){
			coluna[imax]=$i;
			++imax;
		}
	}
	if(maximo<imax)maximo=imax

	}END{

	imeas=(maximo-int(maximo/10))/colunas;
	for(j=0; j<colunas;++ j){
		media[j]=0.0;
		media2[j]=0.0;
		sets[j]=0.0;
	}
	for(i=imeas; i<maximo/colunas; ++i){
		for(j=0; j<colunas;++ j){
			media[j]+=coluna[i*colunas+j];
			media2[j]+=coluna[i*colunas+j]*coluna[i*colunas+j];
			++sets[j];
		}
	}
	for(j=0; j<colunas;++ j){
		media[j]/=sets[j];
		media2[j]/=sets[j];
		if(j==0){
			std=media2[j]-media[j]*media[j];
			std=sqrt(std);
		}
	}
	for(j=0; j<colunas; ++j){
		if(j==0){
			printf("%f ",media[1]);
		}else if(j==1){
			printf("%f %f ",media[0],std);
		}else{
			printf("%f ",media[j]);
		}
	}
	printf("\n");
}
