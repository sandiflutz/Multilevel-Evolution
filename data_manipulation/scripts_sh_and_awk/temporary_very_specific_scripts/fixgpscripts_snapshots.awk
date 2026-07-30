{
	if(FNR==1){
		printf("unset tics\n");
		printf("set size square\n");
		printf("unset border\n");
		printf("set lmargin 0\n");
		printf("set rmargin 0\n");
		printf("set bmargin 0\n");
		printf("set tmargin 0\n");
		printf("unset border\n");
	}
	for(i=1; i<FNR; ++i){
		printf("%s",$i);
	}
	printf("\n");
}
