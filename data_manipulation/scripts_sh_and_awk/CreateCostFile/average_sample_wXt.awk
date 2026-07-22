{
 	if(FNR==1){
              imax=0;
	      ++files;
	}
	if($1 !~ /#/){
		if (NF>colunas) colunas = NF;
    		for(i=1; i <=  NF; ++i){
	 		coluna[imax]  += $i;
         		++sets[imax]; 
	 		++imax;
		} 
    	}
# If it's a comment, copy it but only once:
    else if (files==1){
            printf("#%s ",paramname);
	    for(i=2; i<=NF; ++i){
		    printf("%s ",$1);
	    }
	    printf("%d:sample %d:time\n",colunas,colunas+1);
    }
# The largest file has maximo/NF lines:
     if (imax>maximo) maximo = imax;

}
END {
     for(i=0; i< maximo/colunas; ++i){
	  printf("%f ",param);
          for(j=1; j < colunas; ++j){
		media = coluna[i*colunas+j]/sets[i*colunas+j];
		printf("%.8f ",media);
	  }
	  printf("%d %f\n",sets[i*colunas],coluna[i*colunas]/sets[i*colunas]);
     }
}


