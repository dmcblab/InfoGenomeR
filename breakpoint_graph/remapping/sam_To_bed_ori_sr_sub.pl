my $fasta=$ARGV[0];

use Bio::DB::Fasta;
use warnings;
my $CIGAR=$ARGV[1];
my $read_length=$ARGV[2];
my $edge_expand_length = $ARGV[3];
my $insert_min=210;
my $insert_max=900;
my $max_mismatch=3;
my $ref_seq = Bio::DB::Fasta->new($fasta);
my $edge_seq = Bio::DB::Fasta->new("edge_filling.fa");
$refm_file="NPE_mapping_to_ref_sorted_paired_primary.sam";
$tarm_file="edge_filling_sorted_paired_primary.sam";
my $phred=33;
#$file="test1.txt";
open(my $refm, '<', $refm_file) or die "no file";
open(my $tarm, '<', $tarm_file) or die "no file";
#open(my $unsatisfied_edge, '<', "unsatisfied.list.no_telomere_ends") or die "no file";
my $refm_line;
my $tarm_line;
my @tarm_fields =();
my @read1=();
my @read2=();
my $read1_seq;
my $read2_seq;

my @total_read=();

while($tarm_line = <$tarm>){
	@total_read=();
	@tarm_fields =();
	@read1=();
	@read2=();
	$read1_seq="";
	$read2_seq="";
	edge_pairing($tarm_line, "tarm");
	edge_pairing_print("tarm");

        @tarm_fields =();
        @read1=();
        @read2=();
	$read1_seq="";
	$read2_seq="";
	$refm_line = <$refm>;
	edge_pairing($refm_line, "refm");
	edge_pairing_print("refm");
        if($#total_read != -1){
                my $total_read_i;
                my $total_read_j;
		my $sum_of_p = 0 ;
                for($total_read_i = 0; $total_read_i <= $#total_read; $total_read_i++){
			$sum_of_p = $sum_of_p + 10**-($total_read[$total_read_i][9] + $total_read[$total_read_i][10]);
		}
		
                for($total_read_i = 0; $total_read_i <= $#total_read; $total_read_i++){
                        my $mapping_p = (10**-($total_read[$total_read_i][9] + $total_read[$total_read_i][10]))/$sum_of_p;;
                        for($total_read_j = 0; $total_read_j <= 10; $total_read_j++){
                                print "$total_read[$total_read_i][$total_read_j]\t";
                        }
			print $mapping_p;
                        print "\n";
                }
        }

}
sub edge_pairing_print{
        if($#read1 !=-1 && $#read2 !=-1){
        for(my $j =0; $j <=$#read1; $j++){
                for(my $k=0; $k<=$#read2; $k++){
                        $pos1=$read1[$j][1];
                        $pos2=$read2[$k][1];
			my $first_pos;
			my $last_pos;
			if(abs ($pos1) < abs ($pos2)){
				$first_pos = $pos1;
				$last_pos = $pos2;
			}else{
				$first_pos = $pos2;
				$last_pos = $pos1;
			}
                        $pos1 =~ s/^\+|^\-//;
                        $pos2 =~ s/^\+|^\-//;
			my $insert_size = abs ($pos1-$pos2) + 100;
			my $chr_same_test=0;
			if($read1[$j][0] eq $read2[$k][0]){
				$chr_same_test = 1;
			}
                        if($_[0] eq "tarm"){
                                $first_pos = 1;
                                $last_pos = -1;
				$insert_size = ($insert_min + $insert_max) /2;
				$chr_same_test =1;
                        }

			if($chr_same_test && $insert_size > ($insert_min) && $insert_size < ($insert_max) && $first_pos > 0 && $last_pos < 0 ){
				if(substr($read1[$j][1],0,1) eq "+"){
					if($_[0] eq "tarm"){$temp1_seq = $edge_seq->seq($read1[$j][0], $pos1=>$pos1+$read_length-1);}else{$temp1_seq = $ref_seq->seq($read1[$j][0], $pos1=>$pos1+$read_length-1);}
				}else{
					if($_[0] eq "tarm"){$temp1_seq = $edge_seq->seq($read1[$j][0], $pos1=>$pos1+$read_length-1);}else{$temp1_seq = $ref_seq->seq($read1[$j][0], $pos1=>$pos1+$read_length-1);}
					$temp1_seq =~ tr/ACGTacgt/TGCAtgca/;
					$temp1_seq = scalar reverse $temp1_seq;
				}

				if(substr($read2[$k][1],0,1) eq "+"){
					if($_[0] eq "tarm"){$temp2_seq = $edge_seq->seq($read2[$k][0], $pos2=>$pos2+$read_length-1);}else{$temp2_seq = $ref_seq->seq($read2[$k][0], $pos2=>$pos2+$read_length-1);}
				}else{
					if($_[0] eq "tarm"){$temp2_seq = $edge_seq->seq($read2[$k][0], $pos2=>$pos2+$read_length-1);}else{$temp2_seq = $ref_seq->seq($read2[$k][0], $pos2=>$pos2+$read_length-1);}
					$temp2_seq =~ tr/ACGTacgt/TGCAtgca/;
					$temp2_seq = scalar reverse $temp2_seq;
				}
				$temp1_seq = uc $temp1_seq;
				$temp2_seq = uc $temp2_seq;

				#print "$temp1_seq\n$read1_seq\n";
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
				$refm_not_exist_in_tarm=1;
######################
#				if($_[0] eq "refm"){
#					open(my $unsatisfied_edge, '<', "unsatisfied.list.no_telomere_ends") or die "no file";
#					my $tarm_exist=0;
#					while($ue_line = <$unsatisfied_edge>){
# 				               chomp $ue_line;
#				               my @ue_fields = split("\t", $ue_line);
#						if($read1[$j][0] eq "chr".$ue_fields[0] && abs (abs ($read1[$j][1]) - $ue_fields[1]) < $edge_expand_length){
#							$tarm_exist=$tarm_exist+1;	
#						}
#                                                if($read2[$k][0] eq "chr".$ue_fields[0] && abs (abs ($read2[$k][1]) - $ue_fields[1]) < $edge_expand_length){
#                                                        $tarm_exist=$tarm_exist+1;
#                                                }
#					}
#					if($tarm_exist == 2){
#						$refm_not_exist_in_tarm=0;
#                    			        print "why\t$read1[$j][0]\t$read1[$j][1]\t$read2[$k][0]\t$read2[$k][1]\t$tarm_fields[0]\t$read1[$j][2]\t$read2[$k][2]\t$nm_count1\t$nm_count2\t$nm_qual_sum1\t$nm_qual_sum2\n";
#                
#					}
#				}
########################################################### WE DONT NEED TO TEST BECAUSE WE ONLY SELECTED CORRECTLY MAPPED READS AT REFM #############

				if($refm_not_exist_in_tarm){
					$total_read[$#total_read+1][0]=$read1[$j][0];
	                                $total_read[$#total_read][1]=$read1[$j][1];
	                                $total_read[$#total_read][2]=$read2[$k][0];
	                                $total_read[$#total_read][3]=$read2[$k][1];
	                                $total_read[$#total_read][4]=$tarm_fields[0];
	                                $total_read[$#total_read][5]=$read1[$j][2];
	                                $total_read[$#total_read][6]=$read2[$k][2];
	                                $total_read[$#total_read][7]=$nm_count1;
	                                $total_read[$#total_read][8]=$nm_count2;
	                                $total_read[$#total_read][9]=$nm_qual_sum1;
	                                $total_read[$#total_read][10]=$nm_qual_sum2;
				}
		#               print "$read1[$j][0]\t$read1[$j][1]\t$read2[$k][0]\t$read2[$k][1]\t$tarm_fields[0]\t$read1[$j][2]\t$read2[$k][2]\t$temp1_seq\t$read1_seq\t$temp2_seq\t$read2_seq\n";
		#		print "$read1[$j][0]\t$read1[$j][1]\t$read2[$k][0]\t$read2[$k][1]\t$tarm_fields[0]\t$read1[$j][2]\t$read2[$k][2]\t$nm_count1\t$nm_count2\t$nm_qual_sum1\t$nm_qual_sum2\n";
			}
	      }
	}
        }
}

sub edge_pairing{
        my $sub_line=$_[0];
	for(my $i=0; $i<=1; $i++){
		my @read=();
		my $read_i=0;
		chomp $sub_line;
		@tarm_fields = split("\t", $sub_line);
                my @nm_t= split("NM:i:", $sub_line);
		if($#nm_t > 0){
			my @nm= split("\t", $nm_t[1]);
			if($tarm_fields[5] eq $CIGAR && $nm[0] <=$max_mismatch){
				if((16 & $tarm_fields[1]) ne 16){
					$read[$read_i][1]="+$tarm_fields[3]";
				}else{
					$read[$read_i][1]="-$tarm_fields[3]";
				}
				$read[$read_i][0]=$tarm_fields[2];
				$read[$read_i][2]=$tarm_fields[5];
				$read_i=$read_i+1;

			}
			@XA_t = split("XA:Z:", $sub_line);
			$XA_n = $#XA_t+1;
			if($XA_n > 1){
				@XA = split("\t", $XA_t[1]);
				@xa_fields = split(";", $XA[0]);
				for($xa_i=0;$xa_i <= $#xa_fields;$xa_i++){
					@xa_value = split(",", $xa_fields[$xa_i]);
					if($xa_value[2] eq $CIGAR && $xa_value[3] <=$max_mismatch){
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
			if($_[1] eq "tarm"){
				$tarm_line = <$tarm>;
				$sub_line = $tarm_line;
			}else{
			        $refm_line = <$refm>;
				$sub_line = $refm_line;
			}
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
}	
