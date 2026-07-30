{
	if(INVBIN!=0){
		invbin=INVBIN;
	}else{
		invbin=0.01;
	}
	INVBIN
	if($2==invbin){
		print $1,$3
	}
}
