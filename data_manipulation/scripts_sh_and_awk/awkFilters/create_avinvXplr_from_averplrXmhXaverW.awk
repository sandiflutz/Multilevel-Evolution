{
	mh=MH/1.0;
	if($1==mh){
		for(i=2; i<=NF; ++i){
			printf("%f ",$i);
		}
		printf("\n");
	}
	
}
