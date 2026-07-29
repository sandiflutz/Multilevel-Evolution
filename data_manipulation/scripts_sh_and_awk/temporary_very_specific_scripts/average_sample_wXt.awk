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
    else if (files==1) print $0
# The largest file has maximo/NF lines:
     if (imax>maximo) maximo = imax;

}
END {
     for(i=0; i< maximo/colunas; ++i){
          for(j=0; j < colunas; ++j){
		  media = coluna[i*colunas+j]/sets[i*colunas+j];
		 # if((j==1)||(j==11)){
		  if((j==meascol)||(j==paramcol)){
			printf("%.8f ",media);
		  }
	  }
	  printf("%d\n",sets[i*colunas]);
     }
}


