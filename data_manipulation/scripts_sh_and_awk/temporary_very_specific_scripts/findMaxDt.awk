{
	if(FNR==1){
		maxdt=0.0;
	} 
	if(($2>maxdt)&&($3>=0.9))maxdt=$2;

}END{
	printf("\n%f\n",maxdt);
}
