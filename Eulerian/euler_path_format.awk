{
	if(NF>1){
		for(i=1; i<=NF;i++){
			j=i+1; 
			if(i==1 ){
				printf "%d  ", $i;
				printf "%d-%d  ", $i, $j;
			}else if(i==NF){
				printf "%d", $i;
			}else{
				printf "%d-%d  ", $i, $j;
			}
		}
	}; 
	printf "\n";
}
