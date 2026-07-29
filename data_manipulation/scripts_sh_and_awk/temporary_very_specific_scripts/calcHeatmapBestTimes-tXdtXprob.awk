BEGIN{
        idt=0;
	idtmax=0;
	wbtr=0.9;
}{
	if ($1 !~ /#/){
		time[idt]=$1;
		dt[idt]=$2;
		wcli[idt]=$3;
		wclf[idt]=$4;
		averwi[idt]=$5;
		averwf[idt]=$6;
		clsizei[idt]=$7;
		clsizef[idt]=$8;

		++idt;
		if(idt>idtmax)idtmax=idt;
	}

	}END{

	printf("#1:time[i+1] 2:dt[i+1]-dt[i] 3:wcli[i+1]-wclf[i] 4:wcli[i+1] 5:averwi[i+1] 6:averwf[i+1] 7:clsizef[i+1]-clsizei[i+1]\n",time[i+1],dt[i+1]-dt[i],wcli[i+1]-wclf[i],wcli[i+1],averwi[i+1],averf[i+1],clsizef[i+1]-clsizei[i+1]);
	for(i=0; i<idtmax-1; ++i){
		printf("%f %f %f %f %f %f %d\n",time[i+1],dt[i+1]-dt[i],wcli[i+1]-wclf[i],wcli[i+1],averwi[i+1],averf[i+1],clsizef[i+1]-clsizei[i+1]);
	}
}
