my $library = $ARGV[0];
my $mean_phy_depth = $ARGV[1];
my $search_length = $ARGV[2];
#`perl edge_filling.pl > edge_filling.fa`;
#`bwa-0.7.15 index edge_filling.fa`;
#`bwa-0.7.15 mem -t 20 -h 20 edge_filling.fa ../../NPE.fq1 ../../NPE.fq2 | samtools view -bS > edge_filling.bam`;
#`samtools sort -n edge_filling.bam > edge_filling_sorted.bam`;
#`samtools view -f 1 -F 3840 edge_filling_sorted.bam > edge_filling_sorted_paired_primary.sam`;
#`bwa-0.7.15 mem -t 20 -h 20 /DATA1/AITL/hg19/ucsc.hg19.fasta ../../NPE.fq1 ../../NPE.fq2 | samtools view -bS > NPE_mapping_to_ref.bam`;
#`samtools sort -n NPE_mapping_to_ref.bam > NPE_mapping_to_ref_sorted.bam`;
#`samtools view -f 1 -F 3840 NPE_mapping_to_ref_sorted.bam > NPE_mapping_to_ref_sorted_paired_primary.sam`;
#`mkdir edge_filling_reads`;
#`perl  sam_To_bed_ori_sr_sub.pl > edge_filling_sorted_paired_primary.sam.info`;
`mkdir -p edge_filling_reads_recal`;

my @unsatisfied;


open($unsatisfied_fh, '<', 'unsatisfied.list.no_telomere_ends');

my $unsatisfied_i = 0;
my @line;
while(my $line = <$unsatisfied_fh>){
        @line = ();
        chomp $line;
        next if $line =~ /^\#/ || $line =~ /^\s*$/ || $line =~ /^\+/;
        @line = split(/\t+/, $line);
	$unsatisfied[$unsatisfied_i][0]=$line[0];
        $unsatisfied[$unsatisfied_i][1]=$line[1];
        $unsatisfied[$unsatisfied_i][2]=$line[2];
        $unsatisfied[$unsatisfied_i][3]=$line[3];
        $unsatisfied[$unsatisfied_i][4]=$line[4];
	$unsatisfied_i=$unsatisfied_i+1;
}
open($read_info, '<', 'edge_filling_sorted_paired_primary.sam.info.recal');

my $read_out;
while(my $line = <$read_info>){
        @line = ();
        chomp $line;
        next if $line =~ /^\#/ || $line =~ /^\s*$/ || $line =~ /^\+/;
        @line = split(/\t+/, $line);
         my @sign1=split //, $line[1];
         my @sign2=split //, $line[3];
        my @node1=split /\_/, $line[0];
        my @node2=split /\_/, $line[2];

        my $read_out;
        if($node1[1]<$node2[1]){
                open($read_out, '>>', "./edge_filling_reads_recal/$line[0].$sign1[0].$line[2].$sign2[0]");
        }else{
                open($read_out, '>>', "./edge_filling_reads_recal/$line[2].$sign2[0].$line[0].$sign1[0]");
        }

        $sign1=join("",@sign1);
        $sign2=join("",@sign2);

        if($node1[1]<$node2[1]){
                print $read_out "$line[0]\t$sign1\t$line[2]\t$sign2\t$line[4]\t$line[7]\t$line[8]\t$line[9]\t$line[10]\t$line[11]\n";
        }else{
                print $read_out "$line[2]\t$sign2\t$line[0]\t$sign1\t$line[4]\t$line[8]\t$line[7]\t$line[10]\t$line[9]\t$line[11]\n";
        }
}

for(my $unsatisfied_i = 0; $unsatisfied_i <$#unsatisfied; $unsatisfied_i++){
	for(my $unsatisfied_j = $unsatisfied_i+1; $unsatisfied_j <=$#unsatisfied; $unsatisfied_j++){
		my $line1_index=$unsatisfied_i+1;
		my $ori1_pre=$unsatisfied[$unsatisfied_i][0] . '_' . $line1_index;
		my $ori1;
		if($unsatisfied[$unsatisfied_i][2] % 2 == 0){
			$ori1="+"
		}else{
			$ori1="-"
		}
                my $line2_index=$unsatisfied_j+1;
                my $ori2_pre=$unsatisfied[$unsatisfied_j][0] . '_' . $line2_index;
                my $ori2;
                if($unsatisfied[$unsatisfied_j][2] % 2 == 0){
                        $ori2="+"
                }else{
                        $ori2="-"
                }
	#       print "\n${ori1_pre}\t$unsatisfied[$unsatisfied_i][2]\t${ori2_pre}\t$unsatisfied[$unsatisfied_j][2]\n";

		 $read_in= "./edge_filling_reads_recal/${ori1_pre}.$ori1.${ori2_pre}.$ori2";
		my @stat = stat $read_in;
		
		if( $stat[7] != 0  && $unsatisfied[$unsatisfied_i][3] eq $unsatisfied[$unsatisfied_j][3]){
			$c=`Rscript $library/remapping/edge_filling_max_cal.R $read_in $unsatisfied[$unsatisfied_i][3] $mean_phy_depth $search_length`;
			@c = split(/\s+/, $c);
			$c = $c[1];
		}else{
			$c = 0;
		}
               print "${ori1_pre}\t$unsatisfied[$unsatisfied_i][2]\t${ori2_pre}\t$unsatisfied[$unsatisfied_j][2]\t$c\n";


	}
}


