#my $final_index_of_chr=$ARGV[0];
my $IND=$ARGV[0];
my $final_index_of_chr=`ls SVs.results.* | awk 'BEGIN{max=0}{split(\$0,f,\".\"); if(max<f[3]){max=f[3]}}END{print max}'`;
#my $final_index_of_chr=`SVs.results.* | awk 'BEGIN{max=0}{split($0,f,".results."); if(max<f[2]){max=f[2]}}END{print max}'`;
use POSIX;
use Bio::DB::Fasta;
use Switch;
my @chr_start;
my @chr_end;
my @chr_cumul;
my @chr_hap;
my @chr_number;


for(my $chr_index = 0; $chr_index <=$final_index_of_chr; $chr_index++){
	open(my $results, '<', "SVs.results.${chr_index}");
	my $results_i=0;
        while(my $line = <$results>){
                chomp $line;
                my @line = split(/\t+/, $line);
                $chr_start[$chr_index][$results_i]=$line[0];
                $chr_end[$chr_index][$results_i]=$line[1];
                $chr_cumul[$chr_index][$results_i]=$line[2];
                $chr_hap[$chr_index][$results_i]=$line[3];
                $chr_number[$chr_index][$results_i]=$line[4];
		$results_i++;
	}
}
#open(my $germ_SVs, '<', "/DASstorage6/leeyh/1000G_SNP_indel_SV_integrated/${IND}.germline_initial_SVs");
open(my $germ_SVs, '<', "${IND}.germline_initial_SVs");

while(my $line = <$germ_SVs>){
	chomp $line;
        my @line = split(/\t+/, $line);
	my @test_coor=();
	if($line[0] eq "3to5"){
		$test_coor[$#test_coor+1] = $line[8];
	}elsif($line[0] eq "5to3"){
		for($cn=1; $cn<=$line[7]; $cn++){
		$test_coor[$#test_coor+1] = $line[8]+($line[6]-$line[3])*$cn;
		}
	}
       # print "$#test_coor\t";
	my $changed_cn = old_coor_exist_test($line[1], $line[2],@test_coor);
	if($changed_cn != 0){
		print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$changed_cn\t$line[8]\n";
	}	
}


sub old_coor_exist_test{
        my ($old_fa_hap, $old_fa_chr, @old_fa_coor)=(@_);
	my @where_is;
#	print "$old_fa_hap\t$old_fa_chr\t";
#	print "$old_fa_coor[0]\n";
	for(my $old_fa_coor_i=0;$old_fa_coor_i<=$#old_fa_coor;$old_fa_coor_i++){
		for(my $new_fa_i=0; $new_fa_i<=$#chr_start; $new_fa_i++){
			for(my $new_fa_j=0; $new_fa_j<=$#{@chr_cumul[$new_fa_i]}; $new_fa_j++){
				if($old_fa_hap == $chr_hap[$new_fa_i][$new_fa_j] && $old_fa_chr == $chr_number[$new_fa_i][$new_fa_j] 
					&& $old_fa_coor[$old_fa_coor_i] >= (abs $chr_start[$new_fa_i][$new_fa_j] >= abs $chr_end[$new_fa_i][$new_fa_j] ? abs $chr_end[$new_fa_i][$new_fa_j] : abs $chr_start[$new_fa_i][$new_fa_j]) && $old_fa_coor[$old_fa_coor_i] <= (abs $chr_start[$new_fa_i][$new_fa_j] >= abs $chr_end[$new_fa_i][$new_fa_j] ? abs $chr_start[$new_fa_i][$new_fa_j] : abs $chr_end[$new_fa_i][$new_fa_j])){
					$where_is[$#where_is+1]=$new_fa_i;
				}
			}
		}
	}
	return ($#where_is+1);
}
