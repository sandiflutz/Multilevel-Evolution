{
	if($1 ~/#/){
		printf("#1:t 2:<w> 3:wb 4:clsizeb 5:<wcl⁻> 6:stdwcl⁻ 7:maxclsize 8:wl⁻ 9:<clsize⁻> 10:stdclsize⁻ 11:nh/N 12:cost 13:mb 14:mh 15:rmig 16:rh\n");
	}else{
		printf("%f %f %f %d %f %f %d %f %f %f %f %f %f %f %d %f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,1,$15);
	}
}
