use Bio::DB::Fasta;
use warnings;
my $ref_seq = Bio::DB::Fasta->new("/DATA1/AITL/hg19/ucsc.hg19.fasta");
my $edge_seq = Bio::DB::Fasta->new("edge_filling.fa");
$refm_file="NPE_mapping_to_ref_sorted_paired_primary.sam";
$tarm_file="edge_filling_sorted_paired_primary.sam";
my $phred=33;
#$file="test1.txt";
open(my $refm, '<', $refm_file) or die "no file";
open(my $tarm, '<', $tarm_file) or die "no file";
my $refm_line;
my $tarm_line;
while($tarm_line = <$tarm>){
	my @tarm_fields =();
	my @read1=();
	my @read2=();
	my $read1_seq;
	my $read2_seq;
	for(my $i=0; $i<=1; $i++){
		my @read=();
		my $read_i=0;
		chomp $tarm_line;
		@tarm_fields = split("\t", $tarm_line);
                my @nm_t= split("NM:i:", $tarm_line);
		if($#nm_t > 0){
			my @nm= split("\t", $nm_t[1]);
			if($tarm_fields[5] eq "101M" && $nm[0] <=3){
				if((16 & $tarm_fields[1]) ne 16){
					$read[$read_i][1]="+$tarm_fields[3]";
				}else{
					$read[$read_i][1]="-$tarm_fields[3]";
				}
				$read[$read_i][0]=$tarm_fields[2];
				$read[$read_i][2]=$tarm_fields[5];
				$read_i=$read_i+1;

			}
			@XA_t = split("XA:Z:", $tarm_line);
			$XA_n = $#XA_t+1;
			if($XA_n > 1){
				@XA = split("\t", $XA_t[1]);
				@xa_fields = split(";", $XA[0]);
				for($xa_i=0;$xa_i <= $#xa_fields;$xa_i++){
					@xa_value = split(",", $xa_fields[$xa_i]);
					if($xa_value[2] eq "101M" && $xa_value[3] <=3){
						$read[$read_i][0]=$xa_value[0];
						$read[$read_i][1]=$xa_value[1];
						$read[$read_i][2]=$xa_value[2];
						$read_i=$read_i+1;
					}
				}
			}
		}
		if($i == 0){
			if((16 & $tarm_fields[1]) ne 16){
	                        $read1_seq=$tarm_fields[9];
				$read1_qual=$tarm_fields[10];
			}else{	
				$read1_seq = $tarm_fields[9];
                                $read1_qual=$tarm_fields[10];
				$read1_seq =~ tr/ACGTacgt/TGCAtgca/;
				$read1_seq = scalar reverse $read1_seq;
                                $read1_qual= scalar reverse $read1_qual;


			}
			@read1=@read;
			$tarm_line = <$tarm>;
		}else{
                        if((16 & $tarm_fields[1]) ne 16){
                                $read2_seq=$tarm_fields[9];
                                $read2_qual=$tarm_fields[10];
                        }else{
                                $read2_seq = $tarm_fields[9];
                                $read2_qual=$tarm_fields[10];
                                $read2_seq =~ tr/ACGTacgt/TGCAtgca/;
                                $read2_seq = scalar reverse $read2_seq;
                                $read2_qual= scalar reverse $read2_qual;
                        }
			@read2=@read;
		}
	}
	if($#read1 !=-1 && $#read2 !=-1){
	for(my $j =0; $j <=$#read1; $j++){
		for(my $k=0; $k<=$#read2; $k++){
			$pos1=$read1[$j][1];
			$pos1 =~ s/^\+|^\-//;
			$pos2=$read2[$k][1];
			$pos2 =~ s/^\+|^\-//;
			if(substr($read1[$j][1],0,1) eq "+"){
				$temp1_seq = $edge_seq->seq($read1[$j][0], $pos1=>$pos1+100);
			}else{
				$temp1_seq =  $edge_seq->seq($read1[$j][0], $pos1=>$pos1+100);
				$temp1_seq =~ tr/ACGTacgt/TGCAtgca/;
				$temp1_seq = scalar reverse $temp1_seq;
			}
                        if(substr($read2[$k][1],0,1) eq "+"){
                        	$temp2_seq = $edge_seq->seq($read2[$k][0], $pos2=>$pos2+100);
			}else{
                                $temp2_seq = $edge_seq->seq($read2[$k][0], $pos2=>$pos2+100);
				$temp2_seq =~ tr/ACGTacgt/TGCAtgca/;
				$temp2_seq = scalar reverse $temp2_seq;
			}

			my $mask = $temp1_seq ^ $read1_seq;
			my $nm_count1=0;
			my $nm_qual_sum1=0;
			while ($mask =~ /[^\0]/g) {
				$nm_count1=$nm_count1+1;
				$nm_qual_sum1=$nm_qual_sum1+ord(substr($read1_qual,$-[0],1))-$phred;
			}

                        $mask = $temp2_seq ^ $read2_seq;
                        my $nm_count2=0;
                        my $nm_qual_sum2=0;
                        while ($mask =~ /[^\0]/g) {
                                $nm_count2=$nm_count2+1;
                                $nm_qual_sum2=$nm_qual_sum2+ord(substr($read2_qual,$-[0],1))-$phred;
                        }



	#		print "$read1[$j][0]\t$read1[$j][1]\t$read2[$k][0]\t$read2[$k][1]\t$tarm_fields[0]\t$read1[$j][2]\t$read2[$k][2]\t$temp1_seq\t$read1_seq\t$temp2_seq\t$read2_seq\n";
                        print "$read1[$j][0]\t$read1[$j][1]\t$read2[$k][0]\t$read2[$k][1]\t$tarm_fields[0]\t$read1[$j][2]\t$read2[$k][2]\t$nm_count1\t$nm_count2\t$nm_qual_sum1\t$nm_qual_sum2\n";

		}
	}
	}
}	
