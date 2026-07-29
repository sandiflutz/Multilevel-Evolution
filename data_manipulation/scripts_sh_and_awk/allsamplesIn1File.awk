{
	if(FNR==1){
		printf("\n");
	}
	if($1!~/#/){
		for(i=1; i<NF; ++i){
			printf("%s ",$i);
		}
		printf("\n");
	}
}
