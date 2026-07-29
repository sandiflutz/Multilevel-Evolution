{
		if (FNR==1) {
			imax=0;
		}

		if ($1 !~ /#/){

	
			if (NF>colunas) colunas = NF;

			for (i=1; i <=  NF; ++i){
				coluna[imax]=$i;
				++imax;
			}
		}else{
			for(i=1; i <=  NF; ++i){
				printf("%s ",$i);
			}
			printf("\n");
		}
		
		if (imax>maximo) maximo = imax;
}END{

	for(j=0; j<colunas; ++j){
		printf("%f ",coluna[j]);
	}
	printf("\n");

	for(i=1; i<maximo/colunas; ++i){
		#first column is written separately
		if(coluna[i*colunas]==coluna[(i-1)*colunas]){
			printf("%f ",coluna[i*colunas]);
		}else{
			printf("\n%f ",coluna[i*colunas]);
		}
		#other columns	
		for(j=1; j<colunas; ++j){
			printf("%f ",coluna[i*colunas+j]);
		}
		printf("\n");

	}
}
