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
	mediadt=0.0;
	mediadt2=0.0;
	if(colunas<maxcol){
		maxcol=colunas;
	}
	for(i=0; i<maximo/colunas; ++i){
		for(j=0; j<maxcol; ++j){
			media[j]+=coluna[i*colunas+j];
			media2[j]+=coluna[i*colunas+j]*coluna[i*colunas+j];
			++sets[j];
		}
		mediadt+=coluna[i*colunas+1]-coluna[i*colunas+2]
		mediadt2+=(coluna[i*colunas+1]-coluna[i*colunas+2])*(coluna[i*colunas+1]-coluna[i*colunas+2]);
	}
	mediadt=mediadt*colunas/maximo;
	mediadt2=mediadt2*colunas/maximo;
	stddt=mediadt2-mediadt;
	stddt=sqrt(stddt);
	printf("%f ",param);
	for(j=0; j<maxcol; ++j){
		media[j]/=sets[j];
		media2[j]/=sets[j];
		std=media2[j]-media[j]*media[j];
		std=sqrt(std);
		printf("%f %f ",media[j],std);
	}
	printf("%f %f\n",mediadt,stddt);
}
