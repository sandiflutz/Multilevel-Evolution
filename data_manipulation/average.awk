# Comments:
# This script works only with ONE SAMPLE per file!!
# ALL files must have the same number of columns, although
# the number of lines may differ.
# The evaluated errors are sample to sample errors.
# Data lines starting with '#' are comments.
# The columns (in the data file) are numbered from 1 to NF
#
# Command line:  awk -f average.awk input > output
#
{
# For each file read (test if we're reading the first record), 
# resets imax (total number of cells) and increases files:
 if (FNR==1) {
              imax=0;
	      ++files;
	     }
# If not a comment, read record: 
 if ($1 !~ /#/)
    {
# The number of columns: should be constant, but since some lines can
# be empty...
    if (NF>colunas) colunas = NF;
    for (i=1; i <=  NF; ++i)
        {
	 coluna[imax]  += $i;
	 coluna2[imax] += $i*$i;
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
# The integer division gives the number of lines (incomplete lines
# are neglected):
     for(i=0; i< maximo/colunas; ++i){
          for(j=0; j < colunas; ++j){
               
		  media = coluna[i*colunas+j]/sets[i*colunas+j];
               
		 # printf("%.8f ",media);
              
		  if(j==1){
		  	erro = coluna2[i*colunas+j]/sets[i*colunas+j]-media*media;
	      
		  	erro=sqrt(erro);
              
		  	printf("%.8f %.8f ",media,erro);
	  	}else{
			printf("%.8f ",media);
		}
	  }
	  printf("%d\n",sets[i*colunas]);
     }
}


