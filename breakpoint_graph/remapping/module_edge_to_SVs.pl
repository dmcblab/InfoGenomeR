use POSIX;
use warnings;
my $search_length=$ARGV[0];
my $min_SV_length=1000;
open(my $edges, '<', "edge_probability_recal_sorted.filtered");

my @queue;
while(my $line=<$edges>){
	chomp $line;
	my @line = split /\t+/, $line;
	
	my @l1= split /\_/, $line[0];
        my @l2= split /\_/, $line[2];
	
	my $coor1= `sed '$l1[1]!d' unsatisfied.list.no_telomere_ends`; chomp $coor1; my @coor1=split /\t+/, $coor1;
        my $coor2= `sed '$l2[1]!d' unsatisfied.list.no_telomere_ends`; chomp $coor2; my @coor2=split /\t+/, $coor2;
 
	my $ori1;
	my $sign1;
	my $ori2;
	my $sign2;

	if($coor1[2] % 2 == 0 ){
		$ori1= "3";
		$sign1="+";
	}else{
		$ori1= "5";
		$sign1="-";
	}
        if($coor2[2] % 2 == 0 ){
                $ori2= "3";
		$sign2="+";
        }else{
                $ori2= "5";
		$sign2="-";
        }

	my $SVtype="NA";
	if($coor1[0] ne $coor2[0]){
		$SVtype="<TRA>";
	}elsif(abs $coor1[1] - $coor2[1] > $min_SV_length){

		if($ori1 eq $ori2){
			$SVtype="<INV>";
		}else{ 
			if($ori1 eq "5"){
				$SVtype="<DUP>";
			}else{
				$SVtype="<DEL>";
			}
		}
	}

	my $offset1 =`cat edge_filling_reads_recal/$line[0].$sign1.$line[2].$sign2 | awk 'BEGIN{n=0;sum=0;}{if(\$10>0.8){n=n+1;sum=sum+substr(\$2,2,length(\$2)-1)}}END{print int(sum/n)-$search_length}'`;
        my $offset2 =`cat edge_filling_reads_recal/$line[0].$sign1.$line[2].$sign2 | awk 'BEGIN{n=0;sum=0;}{if(\$10>0.8){n=n+1;sum=sum+substr(\$4,2,length(\$4)-1)}}END{print int(sum/n)-$search_length}'`;
	
	$coor1[1] = $coor1[1]+$offset1;
	if ($ori1 == "3"){
		$coor1[1] = $coor1[1]+125;
	}else{
		$coor1[1] = $coor1[1]-125;
	};
	$coor2[1] = $coor2[1]+$offset2;
        if ($ori2 == "3"){
                $coor2[1] = $coor2[1]+125;
        }else{
                $coor2[1] = $coor2[1]-125;
        };

#	$coor1[1]=$coor1[1]+$line[1]-$search_length;
#	$coor2[1]=$coor2[1]+$line[3]-$search_length;
	if($SVtype ne "NA" && in_queue($line[0], $line[2])!=1){
		print "$SVtype\t$coor1[0]\t$coor1[1]\t$coor2[0]\t$coor2[1]\t$ori1", "to$ori2\t", "60\t60\t60\t60\t0\t0\t1\t1\n";
	}
        $queue[$#queue+1]=$line[0];
        $queue[$#queue+1]=$line[2];

}




sub in_queue{
	if($#queue==-1){
		return 0;
	}
	my ($e1, $e2)=(@_);
	for(my $i=0; $i<=$#queue; $i++){
		if($e1 eq $queue[$i] || $e2 eq $queue[$i]){
			return 1;
		}
	}
	return 0;
}
