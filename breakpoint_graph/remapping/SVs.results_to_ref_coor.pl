my $IND=$ARGV[0];
my $type=$ARGV[1];

my $total;
if($type eq "total"){
	$total=`ls SVs.results.* | awk 'BEGIN{max=0}{split(\$0,f,\".\"); if(max<f[3]){max=f[3]}}END{print max}'`;;
}else{
        $total=`ls germline_SVs.results.* | awk 'BEGIN{max=0}{split(\$0,f,\".\"); if(max<f[3]){max=f[3]}}END{print max}'`;;

}
#my $G1000_dir="/DASstorage6/leeyh/1000G_SNP_indel_SV_integrated";



for(my $i=0; $i<=$total; $i++){
	my $results;
	my $results_out;
	if($type eq "total"){
		open($results, '<', "SVs.results.${i}");
	        open($results_out, '>', "SVs.results.${i}.ref_coor");
	}else{
                open($results, '<', "germline_SVs.results.${i}");
                open($results_out, '>', "germline_SVs.results.${i}.ref_coor");
	}

	while(my $line=<$results>){
        	chomp $line;
	        my @line = split(/\t+/, $line);
                my $abs_sign=0;
                if($line[0]<0){
                        $line[0]=abs $line[0];
                        $line[1]=abs $line[1];
                        $abs_sign=1;
                }

		for(my $j=0; $j<=1; $j++){
                open($deltafile, '<', "${IND}.g$line[3].fa.$line[4].delta") or die;

#	        open($deltafile, '<', "/DAS_Storage5/leeyh/2017.08.04.backup_from_NGS/BG_job/2017.10.3.simulations_from_NGS/10.g$line[3].fa.$line[4].delta") or die;
 	        my @delta_line;
	        my $delta_line = <$deltafile>;
	        my $delta=0;
	        my $last_found=0;

		        while($delta_line = <$deltafile>){
		                @delta_line = ();
		                chomp $delta_line;
		                @delta_line = split(/\t+/, $delta_line);
	 	               if($line[$j]<=$delta_line[0]){
		                        $line[$j]=$line[$j]+$delta;
		                        $last_found=1;
		                        last;
				}
	 	               $delta = $delta_line[1];
	 	       }close($deltafile);
	 	       if($last_found==0){
		                $line[$j]=$line[$j]+$delta;
		        }
		}

		if($abs_sign==1){
			$line[0]=-$line[0];
			$line[1]=-$line[1];
		}


 	       print $results_out "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\n";

	}


}

