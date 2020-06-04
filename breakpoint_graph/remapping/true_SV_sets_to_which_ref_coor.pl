my $IND=$ARGV[0];
#my $G1000_dir="/DASstorage6/leeyh/1000G_SNP_indel_SV_integrated";

open(my $SVs, '<', "$ARGV[1]");
open(my $SVs_post, '>', "$ARGV[2]");
while($line=<$SVs>){
	chomp $line;
        my @line = split(/\t+/, $line);
	my $deltafile;
        open($deltafile, '<', "${IND}.g$line[1].fa.$line[2].delta") or die;
	my @delta_line;
	my $delta_line = <$deltafile>;
	my $delta=0;
	my $last_found=0;
	while($delta_line = <$deltafile>){
		@delta_line = ();
		chomp $delta_line;
		@delta_line = split(/\t+/, $delta_line);
		if($line[3]<=$delta_line[0]){
			$line[3]=$line[3]+$delta;
			$last_found=1;
			last;
		}
		$delta = $delta_line[1];
	}close($deltafile);
	if($last_found==0){
		$line[3]=$line[3]+$delta;
	}

        open($deltafile, '<', "${IND}.g$line[4].fa.$line[5].delta") or die;
	$delta_line = <$deltafile>;
	$delta=0;
	$last_found=0;
	while($delta_line = <$deltafile>){
		@delta_line = ();
		chomp $delta_line;
		@delta_line = split(/\t+/, $delta_line);
		if($line[6]<=$delta_line[0]){
			$line[6]=$line[6]+$delta;
			$last_found=1;
			last;
		}
		$delta = $delta_line[1];
	}close($deltafile);
	if($last_found==0){
		$line[6]=$line[6]+$delta;
	}

	print $SVs_post  "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$line[6]\t$line[7]\n";
}
