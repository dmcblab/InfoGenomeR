my $fasta=$ARGV[0];
my $chr=$ARGV[1];
my $start=$ARGV[2];
my $end=$ARGV[3];
my $reverse_com=$ARGV[4];
use Bio::DB::Fasta;
use warnings;
my $ref_seq = Bio::DB::Fasta->new($fasta);
my $seq= $ref_seq->seq($chr,$start=>$end);

if ( $reverse_com eq "T"){
	$seq =~ tr/ACGTacgt/TGCAtgca/;
        $seq = scalar reverse $seq;
	print "$seq\n";
}else{
	print "$seq\n";
}


