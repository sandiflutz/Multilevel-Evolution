{
	if (MULT>0){
		mult=MULT;
	}else{
		mult=1;
	}
	if (IDCOL!=0){
		idcol=IDCOL;
	}else{
		idcol=1;
	}
	if(FNR==1) {
              ++files;
        }
	if ($1 !~ /#/){
	
		if (MAXCOLS!=0){
			colunas=MAXCOLS;
		}else{
			if (NF>colunas) colunas = NF;
		}
		id=int($idcol*mult);
		coluna[idcol-1,id]=$idcol;
		sets[idcol-1,id]=1;
		if(idmax<id)idmax=id; 
    		for (i=1; i <=  colunas; ++i){
			if(i!=idcol){
         			coluna[i-1,id]  += $i;
         			++sets[i-1,id];
			}
        	}
    	}else if (NR==1){
		for (i=1; i <=  NF; ++i){
			      
			printf("%s ",$i);
		}
		printf("files\n");
	}


}END{

	for(i=0; i<=idmax+1; ++i){
		if(sets[0,i]!=0){
			for(j=0; j<colunas; ++j){
				printf("%f ",coluna[j,i]/sets[j,i]); 
			}
			printf("%d\n",sets[colunas-1,i]);
		}
	}

}
