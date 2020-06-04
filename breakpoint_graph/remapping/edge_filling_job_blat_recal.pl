use warnings;
my $edge_prob_min=0.8;
my $fasta_prefix=$ARGV[0];

`cat edge_probability | sort -k 5,5  > edge_probability_sorted`;
`cat edge_probability_sorted | awk '{if(\$5>0.0000001) print \$0}' > edge_probability_sorted.high`;
open(my $reads_out, '>', "./reads");
open($reads_out, '>>', "./reads");
open(my $edges, '<', 'edge_probability_sorted.high');
while(my $line = <$edges>){
        @line = ();
        chomp $line;
        next if $line =~ /^\#/ || $line =~ /^\s*$/ || $line =~ /^\+/;
        @line = split(/\t+/, $line);
	my $ori1;
	my $ori2;
        if($line[1] % 2 == 0){
        	$ori1="+"
        }else{
        	$ori1="-"
        }
        if($line[3] % 2 == 0){
                $ori2="+"
        }else{
                $ori2="-"
        }
	open(my $reads, '<', "./edge_filling_reads/$line[0].$ori1.$line[2].$ori2");

	while(my $read_line = <$reads>){
		@read_line=();
		chomp $read_line;
		next if $read_line =~ /^\#/ || $read_line =~ /^\s*$/ || $read_line =~ /^\+/;
		@read_line = split(/\t+/, $read_line);
		if($read_line[9] > $edge_prob_min){
			print $reads_out "$read_line[4]\n";
		}
	}

}

`cat edge_filling_sorted_paired_primary.sam | grep -f reads > reads.seq`;

open(my $reads_seq, '<', "./reads.seq");
open(my $reads_seq_arranged, '>', "./reads.seq.arranged");
open(my $reads_seq_arranged_fq, '>', "./reads.seq.arranged.fq");


while(my $seq_line = <$reads_seq>){
        @seq_line = ();
	@seq1 = ();
        chomp $seq_line;
        next if $seq_line =~ /^\#/ || $seq_line =~ /^\s*$/ || $seq_line =~ /^\+/;
        @seq_line = split(/\t+/, $seq_line);
	
	$seq1[0][0]=$seq_line[0];
	$seq1[0][1]=$seq_line[1];
        $seq1[0][2]=$seq_line[9];
        $seq1[0][3]=$seq_line[10];

        if((16 & $seq1[0][1]) eq 16 ){
                $seq1[0][2] =~ tr/ACGTacgt/TGCAtgca/;
                $seq1[0][2] = scalar reverse $seq1[0][2];
        	$seq1[0][3] = scalar reverse $seq1[0][3];
        }

	$seq_line = <$reads_seq>;
        @seq_line = ();
        chomp $seq_line;
        next if $seq_line =~ /^\#/ || $seq_line =~ /^\s*$/ || $seq_line =~ /^\+/;
        @seq_line = split(/\t+/, $seq_line);

        $seq1[1][0]=$seq_line[0];
        $seq1[1][1]=$seq_line[1];
        $seq1[1][2]=$seq_line[9];
	$seq1[1][3]=$seq_line[10];
        if((16 & $seq1[1][1]) eq 16 ){
                $seq1[1][2] =~ tr/ACGTacgt/TGCAtgca/;
                $seq1[1][2] = scalar reverse $seq1[1][2];
                $seq1[1][3] = scalar reverse $seq1[1][3];
        }

	if((64 & $seq1[0][1]) ne 64){
 	       $seq1[2][0]=$seq1[0][0];
	       $seq1[2][1]=$seq1[0][1];
	       $seq1[2][2]=$seq1[0][2];
               $seq1[0][0]=$seq1[1][0];
               $seq1[0][1]=$seq1[1][1];
               $seq1[0][2]=$seq1[1][2];
               $seq1[1][0]=$seq1[2][0];
               $seq1[1][1]=$seq1[2][1];
               $seq1[1][2]=$seq1[2][2];
	}

        $seq1[1][2] =~ tr/ACGTacgt/TGCAtgca/;
        $seq1[1][2] = scalar reverse $seq1[1][2];
        $seq1[1][3] = scalar reverse $seq1[1][3];


	print $reads_seq_arranged ">$seq1[0][0]\n";
	print $reads_seq_arranged "$seq1[0][2]";
        print $reads_seq_arranged "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
	print $reads_seq_arranged "$seq1[1][2]\n";

        print $reads_seq_arranged_fq ">$seq1[0][0]\n";
        print $reads_seq_arranged_fq "$seq1[0][2]";
        print $reads_seq_arranged_fq "$seq1[1][2]\n";

        print $reads_seq_arranged_fq "$seq1[0][3]";
        print $reads_seq_arranged_fq "$seq1[1][3]\n";


}

`blat -minMatch=4 $fasta_prefix.2bit reads.seq.arranged output.psl`;
`cat output.psl | awk '{if(NR>5 && \$16-\$17 < 1000 && \$16-\$17 > -1000 && \$1 > 194) print \$0}' > output.psl.filtered`;
