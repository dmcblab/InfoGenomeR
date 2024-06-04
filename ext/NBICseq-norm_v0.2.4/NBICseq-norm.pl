#!/usr/bin/env perl
use strict;

use FindBin qw($Bin);
my $path = $Bin;

my $brs = "$path/BRS_1bp/BRS";
my $randomsample = "$path/randomSample/random_split";
my $prepPois = "$path/DataPrepare/PrepPois";
my $prepPoisGAM = "$path/DataPrepareGAM/PrepPoisGAM";
my $rnormalize = "$path/R/normalize.R";
my $rrefine = "$path/R/refine.R";
my $rrefineGAM = "$path/R/refineGAM.R";
my $rplot_function = "$path/R/plot_RC_vs_GC.R";
my $compRatio = "$path/R/compRatio.R";
my $markCNV = "$path/R/purity.R";
my $purityEM = "$path/purityEM/purityEM";
my $filterCNV = "$path/Filter/filterCNV";
my $rnormalizeGAM = "$path/R/normalizeGAM.R";
my $testmgcv = "$path/R/test.mgcv.installed.R";


if(!(-e $brs)){die"Cannot find the file $brs\n";}
if(!(-e $randomsample)){die "Cannot find the file $randomsample\n";}
if(!(-e $prepPois)){die "Cannot find the file $prepPois\n";}
if(!(-e $prepPoisGAM)){die "Cannot find the file $prepPois\n";}
if(!(-e $rrefine)){ die "Cannot find the file $rrefine\n";}
if(!(-e $rrefineGAM)) {die "Cannot find the file $rrefineGAM\n";}
if(!(-e $rplot_function)) {die "Cannot find the file $rplot_function\n";}
if(!(-e $compRatio)){die  "Cannot find the file $compRatio\n";}
if(!(-e $markCNV)){die "Cannot find the file $markCNV\n";}
if(!(-e $purityEM)){die "Cannot find the file $purityEM\n";}
if(!(-e $filterCNV)){die "Cannot find the file $filterCNV\n";}
if(!(-e $rnormalizeGAM)) {die "Cannot find the file $rnormalizeGAM\n";}
if(!(-e $rnormalize)) {die "Cannot find the file $rnormalize\n";}
if(!(-e $testmgcv)){die"Cannot find the file $testmgcv\n";}


use Getopt::Long;


my $output;
my $readlen=50;
my $fragment_size=300;
my $tmpdir="$path/tmp/";
my $perc = 0.0002;
my $help;
my $gapfile;
my $gc_bin;
my $bin_size;
my $nomapbin;
my $map_bin;
my $bin_only;
my $figname;
my $figtitle="";
my $uds = 5;
my $noGapInRead;
my $filter;

my $invalid;
$invalid = GetOptions("help"=>\$help, "l=i"=>\$readlen, "s=i"=>\$fragment_size, "tmp=s"=>\$tmpdir, "p=f"=>\$perc,"gc_bin"=>\$gc_bin, "b=i"=>\$bin_size, "NoMapBin"=>\$nomapbin, "bin_only"=>\$bin_only, "fig=s"=>\$figname, "title=s"=>\$figtitle, "uds=i"=>\$uds, "noGapInRead"=>\$noGapInRead, "filter"=>\$filter);

my $size = $#ARGV+1;

if($help|!$invalid|| ($size!=1 && $size!=2)) {
        print "Usage: NBICseq-norm.pl [options] <configFile> <output>\n";
        print "Options:\n";
        print "        --help\n";
	print "        -l=<int>: read length\n";
	print "        -s=<int>: fragment size\n";
	print "        -p=<float>: a subsample percentage: default 0.0002.\n";
	print "        -b=<int>: bin the expected and observed as <int> bp bins; Default 100.\n";
	print "        --gc_bin: if specified, report the GC-content in the bins\n";
	#print "        --uds=<int>: the number of nucleotide extended to upstream or downstream; Default 5\n";
	print "        --NoMapBin: if specified, do NOT bin the reads according to the mappability\n";
	#print "        --noGapInRead: if specified, use all nucleotides in the entire extended read in the model traning and prediction\n";
	#print "        --filter: filter the likely CNV regions before run normalization\n";
	print "        --bin_only: only bin the reads without normalization\n";
	print "        --fig=<string>: plot the read count VS GC figure in the specified file (in pdf format)\n";
	print "        --title=<string>: title of the figure\n";
	print "        --tmp=<string>: the tmp directory; If unspecified, use $path/tmp/\n";
	print "        <output> will store the parameter estimate in the Negative Binomial model\n";
	die("\n");
        }

if(!$bin_only && $size!=2){
	print "Must specify <output> when --bin_only is not specified\n";
	die("\n");
	}

if($gc_bin) {$gc_bin = "--gc_bin";} else {$gc_bin="";}
if($nomapbin){$map_bin = "";} else{$map_bin = "--map_bin";}
if($uds<0){ die("The option --uds must be nonnegative\n"); }
if($noGapInRead){
	$noGapInRead = "--noGapInRead";
	}else{
	$noGapInRead = "";
	}


if($tmpdir !~/\/$/){$tmpdir = $tmpdir."\/";}

if(!(-e $tmpdir)){mkdir $tmpdir or die "Cannot create the dirctory $tmpdir\n";}

if($readlen<0) {die "Read Length should be nonnegative\n";}
if($fragment_size<=0) {die "Fragment size should be positive\n";}

if($perc<=0|| $perc>1){
	die "the subsampleing percentage must be between 0 and 1\n";
	}
if($bin_only){
	if(!$bin_size){$bin_size = 10000;}
	}else{
	if(!$bin_size){$bin_size = 100;}
	}
if($bin_size<=0){die("Bin size must be a positive number\n");}

my $config = $ARGV[0];
my $output;
if(!$bin_only) {
	$output = $ARGV[1];
	open(TEMPOUT,">$output") or die("Failed to open the file $output\n");
	close(TEMPOUT)
	}

if(system("R --slave < $testmgcv")!=0){
	die("Please install the R package \'mgcv\'\n");
	}


use File::Temp qw/tempfile/;

open(CONFIG, $config) or die "No such file or directory: $config\n";


#### first try to see if the configure file is correct

my $i = 1;
while(<CONFIG>){
	chomp;
	my @row = split(/\t/);
	my $rowsize = $#row+1;
	if($i>1&& $rowsize==5){
		my $chromname = $row[0];
		my $fa_file = $row[1];
		my $mapfile = $row[2];
		my $readfile = $row[3];
		my $outfile_chr = $row[4];

		if(!(-e $fa_file)) {die "No such file or directory: $fa_file\n";}
		if(!(-e $mapfile)) {die "No such file or directory: $mapfile\n";}
		if(!(-e $readfile)){die "No such file or directory: $readfile\n";}
		if($outfile_chr =~ /,/){die "Output file name should NOT contain comma\n";}

		open FILE, ">$outfile_chr" or die "$! $outfile_chr\n";
		close(FILE);
		}
	$i = $i+1;
	}

close(CONFIG);

#### required to just bin the reads

if($bin_only){
	open(CONFIG, $config) or die "No such file or directory: $config\n";
	
	my $i = 1;
	my $flag = 1;
	my %readcntfile_chr = {};

	my $rplot_1st_arg = "";
	while(<CONFIG>){
		chomp;
		my @row = split(/\t/);
		my $rowsize = $#row+1;
		if($i>1&& $rowsize==5){
			my $chromname = $row[0];
			my $fa_file = $row[1];
			my $mapfile = $row[2];
			my $readfile = $row[3];
			my $outfile_chr = $row[4];

			my $cmd1;
			my $cmd2;


			if(!(exists $readcntfile_chr{$chromname})){
				my $tmp_readcntfile_chr;
				(undef, $tmp_readcntfile_chr) = tempfile("readcount_".$chromname."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);
				$readcntfile_chr{$chromname} = $tmp_readcntfile_chr;
				$cmd1 = "$brs -o $tmp_readcntfile_chr -m $mapfile $readfile";
				$cmd2 = "$prepPois -i $tmp_readcntfile_chr --gc_bin $map_bin -b $bin_size -o $outfile_chr  $fa_file";
	
				print "$cmd1\n";
				if(system($cmd1)!=0) {die("\n");};
				print "$cmd2\n\n";
				if(system($cmd2)!=0) {die("\n");};
				$rplot_1st_arg = $rplot_1st_arg."$outfile_chr".",";

				system("rm $tmp_readcntfile_chr");
				}
			}
		$i = $i +1;
		}

	if($rplot_1st_arg =~ /,$/) {chop($rplot_1st_arg);}
	if($figname){
		my $cmd = "R --slave --args $rplot_1st_arg $figname \"$figtitle\" < $rplot_function";
		print "$cmd\n";
		if(system($cmd)!=0){die("\n");}
		}

	exit;
	}


##########################################################################################
###################### then run normalization ############################################
##########################################################################################
open(CONFIG, $config) or die "No such file or directory: $config\n";


my $tmpfile_total_reads;
(undef, $tmpfile_total_reads) = tempfile("Summary_ReadCount"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);

my $tmpfile_total_reads_tmp;
(undef, $tmpfile_total_reads_tmp) = tempfile("Summary_ReadCount_tmp"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);



###########################################################################################
####################### first create the 1 bp bins ########################################
###########################################################################################
my $i = 1;
my $flag = 1;
my %readcntfile_chr = {};
while(<CONFIG>){
        chomp;
        my @row = split(/\t/);
        my $rowsize = $#row+1;
        if($i>1&& $rowsize==5){
                my $chromname = $row[0];
                my $fa_file = $row[1];
                my $mapfile = $row[2];
                my $readfile = $row[3];

		my $cmd1;
		if(!(exists $readcntfile_chr{$chromname})){
			my $tmp_readcntfile_chr;
			(undef, $tmp_readcntfile_chr) = tempfile("readcount_".$chromname."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);
			$readcntfile_chr{$chromname} = $tmp_readcntfile_chr;
			if($flag==1) {
				$cmd1 = "$brs -o $tmp_readcntfile_chr -m $mapfile -s $tmpfile_total_reads_tmp $readfile";
				$flag = 0;
				}else{
                                $cmd1 = "$brs -o $tmp_readcntfile_chr -m $mapfile -s $tmpfile_total_reads_tmp $readfile";
				}
			print "$cmd1\n";
			if(system($cmd1)!=0){ die("\n");}
			my $cmd3 = "tail -1 $tmpfile_total_reads_tmp >> $tmpfile_total_reads";
			print "$cmd3\n";
			system($cmd3);
			}else{
			print STDERR "Warning: two $chromname in the configure file $config\nOnly used the first one\n";
			}
                }
        $i = $i+1;
        }
close(CONFIG);


open(FILE_TREAD ,$tmpfile_total_reads) or die "No such file or directory: $tmpfile_total_reads\n";
my $total_read = 0;
my $total_mappable_pos = 0;
while(<FILE_TREAD>){
        chomp;
        my @row = split(/\t/);
        $total_read = $total_read + $row[0];
        $total_mappable_pos = $total_mappable_pos + $row[1];
        }
close(FILE_TREAD);

if(-e $tmpfile_total_reads) {system("rm $tmpfile_total_reads");}
if(-e $tmpfile_total_reads_tmp) {system("rm $tmpfile_total_reads_tmp");}

### estimate the bin size used for the region filtering step and the refinement step
use POSIX;
my $bin_size_1ststep = ceil(1000/($total_read/$total_mappable_pos));
if($bin_size_1ststep < 1000) {$bin_size_1ststep = 1000;}
print "bin size is $bin_size_1ststep\n";

##############################################################################################
######## now bin the reads and try to filter the regions that are likely CNV regions #########
##############################################################################################

my $tmp_filter_binfile;
my $tmp_filter_ratiofile; ## copy ratio file for purityEM
my $tmp_purityEstimateFile; ## purity and ploidy estimate file
my $tmp_cnvRegionFile; ## cnv region file

#### filter the CNV regions if required
if($filter){
	(undef, $tmp_filter_binfile) = tempfile("binfiltering"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);
	open(CONFIG, $config) or die "No such file or directory: $config\n";
	my $i=1;
	my $flag = 1;
	while(<CONFIG>){
	        chomp;
       		my @row = split(/\t/);
        	my $rowsize = $#row+1;
        	if($i>1&& $rowsize==5){
                	my $chromname = $row[0];
                	my $fa_file = $row[1];
                	my $mapfile = $row[2];
                	my $readfile = $row[3];

                	my $cmd;
                	if(exists $readcntfile_chr{$chromname}){
                        	my $tmp_readcntfile_chr = $readcntfile_chr{$chromname};
                        	if($flag==1) {
					$cmd = "$prepPois -i $tmp_readcntfile_chr --chrom $chromname  --gc_bin --map_bin -b $bin_size_1ststep  $fa_file > $tmp_filter_binfile";
                                	$flag = 0;
                                	}else{
					$cmd = "$prepPois -i $tmp_readcntfile_chr --NoHeader --chrom $chromname --gc_bin --map_bin -b $bin_size_1ststep  $fa_file >> $tmp_filter_binfile";
                                	}
                        	print "$cmd\n";
                        	if(system($cmd)!=0) {die("\n");}
                        	}
                	}
        	$i = $i+1;
        	}
	close(CONFIG);

	### calculate the ratio
	(undef, $tmp_filter_ratiofile) = tempfile("binfiltering"."_ratio_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);
	my $cmd = "R --slave --args $tmp_filter_binfile $tmp_filter_ratiofile < $compRatio";
	### after this step, the file $tmp_filter_binfile will have an aditional column recording copy ratio information
	print "$cmd\n";
	if(system($cmd)!=0){die("\n")};

	### estimate the purity and ploidy
	(undef, $tmp_purityEstimateFile) = tempfile("purityEstimate"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);
	my $cmd = "$purityEM --maxComp 15 --nRS 10 --subsample -o $tmp_purityEstimateFile $tmp_filter_ratiofile\n";
	print "$cmd\n";
	if(system($cmd)!=0){die("\n")};

	### get the potential CNV region
	(undef, $tmp_cnvRegionFile) = tempfile("cnvRegion"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);
	my $cmd = "R --slave --args $tmp_filter_binfile $tmp_purityEstimateFile $tmp_cnvRegionFile $figname < $markCNV";
	print "$cmd\n";
	if(system($cmd)!=0){die("\n")};
	}


#### now the file $tmp_filter_binfile should have an aditional column with copy number estimate


#my $cmd = "rm $tmp_filter_binfile";
# "$cmd\n";print
#die("\n");

#################################################################################
################ prepare data for the GAM  ######################################
#################################################################################

open(CONFIG, $config) or die "No such file or directory: $config\n";

my $tmpfile;
(undef, $tmpfile) = tempfile("poismodel"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);

my $i = 1;
my $flag = 1;
while(<CONFIG>){
        chomp;
        my @row = split(/\t/);
        my $rowsize = $#row+1;
        if($i>1&& $rowsize==5){
                my $chromname = $row[0];
                my $fa_file = $row[1];
                my $mapfile = $row[2];
                my $readfile = $row[3];

                my $cmd1;
                if((exists $readcntfile_chr{$chromname})){
			my $tmp_readcntfile_chr = $readcntfile_chr{$chromname};
                        if($flag==1) {
				if($filter){
					$cmd1 = "$filterCNV $tmp_cnvRegionFile $tmp_readcntfile_chr $chromname | $randomsample $perc | $prepPoisGAM -l $readlen -s $fragment_size --uds $uds $noGapInRead  $fa_file > $tmpfile";
					}else{
					$cmd1 = "cat $tmp_readcntfile_chr | $randomsample $perc | $prepPoisGAM -l $readlen -s $fragment_size --uds $uds $noGapInRead  $fa_file > $tmpfile";
					}
                                $flag = 0;
                                }else{
				if($filter){
					$cmd1 = "$filterCNV $tmp_cnvRegionFile $tmp_readcntfile_chr $chromname | $randomsample $perc | $prepPoisGAM --NoHeader -l $readlen -s $fragment_size --uds $uds $noGapInRead $fa_file >> $tmpfile";
					}else{
					$cmd1 = "cat $tmp_readcntfile_chr | $randomsample $perc | $prepPoisGAM --NoHeader -l $readlen -s $fragment_size --uds $uds $noGapInRead $fa_file >> $tmpfile";
					}
                                }
                        print "$cmd1\n";
			if(system($cmd1)!=0){die("\n")};
                        }
                }
        $i = $i+1;
        }
close(CONFIG);


### use GAM to normalize the data
my $gam_estimate_file;
(undef, $gam_estimate_file) = tempfile("gamModel_estimate"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);
#my $cmd = "R --slave --args $tmpfile $gam_estimate_file < $rnormalize";
my $cmd = "R --slave --args $tmpfile $fragment_size $gam_estimate_file < $rnormalizeGAM";
print "$cmd\n";
if(system($cmd)!=0){die("\n")};




################################################################################################
############################# calculate the expected (the 1st step)#############################
################################################################################################
my $tmp_estimates_file;
(undef, $tmp_estimates_file) = tempfile("parameter_estimate_"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);

my $cmd = "cut -f2,3,4,5 $gam_estimate_file > $tmp_estimates_file";
print "$cmd\n";
if(system($cmd)!=0){die("\n")};

open(CONFIG, $config) or die "No such file or directory: $config\n";
my $i = 1;
my $flag=1;
my $expected_cnt_file_tmp;
(undef, $expected_cnt_file_tmp) = tempfile("count_normalized_1stStep_"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);
while(<CONFIG>){
        chomp;
        my @row = split(/\t/);
        my $rowsize = $#row+1;

        if($i>1&& $rowsize==5){
                my $chromname = $row[0];
                my $fa_file = $row[1];
                my $mapfile = $row[2];
                my $readfile = $row[3];
                #my $expected_cnt_file = $row[4];

                if(exists $readcntfile_chr{$chromname}){
                        my $readcnt_file_tmp = $readcntfile_chr{$chromname};
			my $cmd;
			if($flag==1){
				$cmd = "$prepPoisGAM -l $readlen -s $fragment_size -i $readcnt_file_tmp -b $bin_size_1ststep --gc_bin --map_bin --uds $uds $noGapInRead -e $tmp_estimates_file $fa_file > $expected_cnt_file_tmp";
				$flag=0;
				}
			else{
				$cmd = "$prepPoisGAM -l $readlen -s $fragment_size -i $readcnt_file_tmp -b $bin_size_1ststep --gc_bin --NoHeader --map_bin --uds $uds $noGapInRead -e $tmp_estimates_file  $fa_file >> $expected_cnt_file_tmp";
				}
			print "$cmd\n";
			if(system($cmd)!=0){die("\n")};
                        }
                }
	$i = $i+1;
        }

close(CONFIG);



##Calculate the refined parameters
my $parameter_file_refine;
(undef, $parameter_file_refine) = tempfile("refine_parameter_estimate_"."_XXXXXXX",SUFFIX=>".txt",DIR => $tmpdir,OPEN=>0);

my $cmd;
if($filter){
	$cmd = "R --slave --args $expected_cnt_file_tmp $parameter_file_refine $fragment_size $tmp_filter_binfile < $rrefineGAM";
	}else{
	$cmd = "R --slave --args $expected_cnt_file_tmp $parameter_file_refine $fragment_size < $rrefineGAM";
	}
print "$cmd\n";
if(system($cmd)!=0){die("\n")};


## now I can calculate the refined expected values
 

my $par_refine = "--refine $parameter_file_refine";
open(CONFIG, $config) or die "No such file or directory: $config\n";
my $i = 1;
my $flag=1;
while(<CONFIG>){
        chomp;
        my @row = split(/\t/);
        my $rowsize = $#row+1;

        if($i>1&& $rowsize==5){
                my $chromname = $row[0];
                my $fa_file = $row[1];
                my $mapfile = $row[2];
                my $readfile = $row[3];
                my $expected_cnt_file = $row[4];

                if(exists $readcntfile_chr{$chromname}){
                        my $readcnt_file_tmp = $readcntfile_chr{$chromname};
                        my $cmd;
                        if($flag==1){
				$cmd = "$prepPoisGAM -l $readlen -s $fragment_size -i $readcnt_file_tmp -b $bin_size $gc_bin $par_refine $map_bin --uds $uds $noGapInRead -e $tmp_estimates_file $fa_file > $expected_cnt_file";
                                $flag=0;
                                }
                        else{
				 $cmd = "$prepPoisGAM -l $readlen -s $fragment_size -i $readcnt_file_tmp -b $bin_size $gc_bin $par_refine $map_bin --uds $uds $noGapInRead -e $tmp_estimates_file  $fa_file > $expected_cnt_file";
                                }
                        print "$cmd\n";
			if(system($cmd)!=0){die("\n")};
                        }
                }
        $i = $i+1;
        }

close(CONFIG);


##################################################################################
########### write the important result to $output ################################
##################################################################################

open(OUTPUT, ">$output") or die "Failed to open the file: $output\n";
print OUTPUT "\#FragmentLength $fragment_size\n";
print OUTPUT "\#ReadLen $readlen\n";
print OUTPUT "\#NumberOfExtendedBp $uds\n";
### first the parameters from the GAM model
open(INPUT, "<$gam_estimate_file") or die "Failed to open the file: $gam_estimate_file\n";
my @rows_param = <INPUT>;
close(INPUT);
my $num_lines = $#rows_param+1;
print OUTPUT "\#Parameter estimates from the negative binomial model\n";
print OUTPUT "\# $num_lines lines\n";
foreach(@rows_param){
	chomp;
	print OUTPUT "$_\n";
	}

#### then the parameters from the refinery step
open(INPUT, "<$parameter_file_refine") or die "Failed to open the file: $parameter_file_refine\n";
my @rows_refine = <INPUT>;
close(INPUT);
my $num_lines = $#rows_refine+1;
print OUTPUT "\#Parameter estimates from the refine step\n";
print OUTPUT "\# $num_lines lines\n";
foreach(@rows_refine){
	chomp;
	print OUTPUT "$_\n";
	}
### then the parameters from the purityEM estimate
if($filter){
	open(INPUT,"<$tmp_purityEstimateFile") or die "Failed to open the file: $tmp_purityEstimateFile\n";
	my @rows_param = <INPUT>;
	close(INPUT);
	my $num_lines = $#rows_param+1;
	print OUTPUT "\#Parameter estimates from the purity, ploidy estimation\n";
	print OUTPUT "\# $num_lines lines\n";
	print OUTPUT "\# BinSizeUsed $bin_size_1ststep\n";
	foreach(@rows_param){
		chomp;
	        print OUTPUT "$_\n";
        	}
	}
close(OUTPUT);


###################################################################################
################# now remove the tmp files ########################################
###################################################################################

### remove the temp read count file generated for each chromosome
foreach my $key (keys %readcntfile_chr){
	if(-e $readcntfile_chr{$key}) {
		#print "rm $readcntfile_chr{$key}\n";
		if(-e $readcntfile_chr{$key}) {system("rm $readcntfile_chr{$key}");}
		}
	}



if(-e $tmpfile) {system("rm $tmpfile");} ### temp file of the data prepared for poisson model (subsampled data)
if(-e $tmp_estimates_file) {system("rm $tmp_estimates_file");} ### temp file of parameter estimate in Poisson model 
if(-e $parameter_file_refine) {system("rm $parameter_file_refine");} ## temp file of the parameter estimate in the refining step
if(-e $expected_cnt_file_tmp) {system("rm $expected_cnt_file_tmp");} ## the count and the expected count file from the 1st step

if(-e $tmp_filter_binfile) {system("rm $tmp_filter_binfile");} ## remove the bin file generated for filtering the CNV region
if(-e $tmp_filter_ratiofile) {system("rm $tmp_filter_ratiofile");} ## the copy ratio file for th program purityEM
if(-e $tmp_cnvRegionFile) {system("rm $tmp_cnvRegionFile");} ## the CNV region file for filtering the CNV region
if(-e $tmp_purityEstimateFile) {system("rm $tmp_purityEstimateFile");} ## the parameters from the purity estimation
if(-e $gam_estimate_file) {system("rm $gam_estimate_file");} ## the parameter estimate from the negative binomial model
