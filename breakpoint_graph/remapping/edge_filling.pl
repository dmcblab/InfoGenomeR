my $fasta=$ARGV[0];
my $search_length=$ARGV[1];
my $chr_prefix=$ARGV[2];
use Bio::DB::Fasta;
use strict;
use warnings;
my $db = Bio::DB::Fasta->new($fasta);

my $file = "unsatisfied.list.no_telomere_ends";
open(my $data, '<', $file) or die "no file";
my $index=1;
while (my $line = <$data>) {
        chomp $line;
        my @fields = split("\t", $line);
	print ">$fields[0]_$index\n";
	my $start = $fields[1]-$search_length;
	my $end = $fields[1]+$search_length;
	my $seq;
	if($chr_prefix eq 1){
		$seq = $db->seq("chr$fields[0]", $start=>$end);
	}else{
                $seq = $db->seq("$fields[0]", $start=>$end);
	}
	while( my $chunk = substr($seq, 0, 80, "")){
	        print "$chunk\n";
	}
	$index=$index+1;

}

