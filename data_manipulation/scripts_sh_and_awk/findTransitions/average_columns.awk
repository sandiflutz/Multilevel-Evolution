BEGIN{
	imax=0;
	}{
	if($1 !~ /#/){
		if(colunas<NF)colunas=NF;
		for(i=1; i<=NF; ++i){
			coluna[imax]=$i;
			++imax;
		}
		if(maximo<imax)maximo=imax
	}

	}END{

	for(j=0; j<colunas; ++j){
		media[j]=0.0;
		media2[j]=0.0;
		sets[j]=0;
	}
	for(i=0; i<maximo/colunas; ++i){
		for(j=0; j<colunas; ++j){
			media[j]+=coluna[i*colunas+j];
			media2[j]+=coluna[i*colunas+j]*coluna[i*colunas+j];
			++sets[j];
		}
	}
	printf("%f ",param);
	for(j=0; j<colunas; ++j){
		media[j]/=sets[j];
		media2[j]/=sets[j];
		std=media2[j]-media[j]*media[j];
		std=sqrt(std);
		if(j==1){
			printf("%f %f ",media[j],std);
		}else{
			printf("%f ",media[j]);
		}
	}
	printf("\n");
}
