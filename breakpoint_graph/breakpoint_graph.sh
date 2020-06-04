#!/bin/bash
mode="somatic"
SV="SVs"
#tmp_folder="$PWD/tmp/"
#min_ploidy="1.5";
#max_ploidy="5";
###########
library="/home/qlalf1457/InfoGenomer/TOTAL"
lambda_ini="1"
lambda_fi="4"
bicseq="~/NBICseq-seg_v0.7.2/NBICseq-seg.pl"
#fasta_prefix="/DATA1/AITL/hg19/ucsc.hg19"
#mean_phy_depth=5;
#########################
search_length=50000;
CIGAR="100M";
read_length=100;
###########################
if [ $# != 14 ]; then
        echo $0 \<mode\> \<sample\> \<cancer_type\> \<min_ploidy\> \<max_ploidy\> \<NPE_dir\> \<tumor_bin_dir\> \<germ_bin_dir\> \<germ_LocSeq_result\> \<fasta_prefix\> \<haplotype_coverage\> \<bam\> \<normal_bam\> \<chr_prefix\>
        exit 1
fi
mode=$1
sample=$2
cancer_type=$3
min_ploidy=$4
max_ploidy=$5
NPE_dir=`readlink -f $6`
tumor_bin_dir=`readlink -f $7`
germ_bin_dir=`readlink -f $8`
germ_LocSeq_result=`readlink -f $9`
fasta_prefix=`readlink -f ${10}`
haplotype_coverage=${11}
#bam=${12}
#normal_bam=${13}
bam=`readlink -f ${12}`
normal_bam=`readlink -f ${13}`
chr_prefix=${14}


bin="$tumor_bin_dir/"
normal_bin="$germ_bin_dir/"
#chr_prefix="1"
i=1;
iter=1;

if [ $mode == "germline" ];then
	lambda_ini=1;
#	bin="$working_dir/bicseq_norm/"
	Rscript $library/SV_break_adjustment.R
	$library/SV_truncated_wgs.sh SVs $bin
	Rscript $library/SV_break_adjust.R
	mkdir tmp;
	tmp_folder="$PWD/tmp/"
        Rscript $library/SV_local_CN_segment_wgs_raw.R $SV $i $bin $tmp_folder $bicseq $lambda_ini $bam $chr_prefix;
	rm -rf tmp;
	exit 1
fi



if [ $mode == "somatic" ];then
	mkdir bicseq_norm
	mkdir bicseq_norm_germ

	for chr in {1..23}; do
		Rscript $library/germ_cnv_filter_v2.R $germ_LocSeq_result $normal_bin $bin $chr
	done
	
	bin="$PWD/bicseq_norm/"
	normal_bin="$PWD/bicseq_norm_germ/"
fi


Rscript $library/SV_break_adjustment.R
$library/SV_truncated_wgs.sh SVs $bin


ABSOLUTE_STABLE="F";

for iter in {1..50}
do
	mkdir iter$iter;
	cd iter$iter;
	cp  ../SVs ./
	mkdir tmp;

	abs_line=50000;
        while [ $abs_line -ge 50000 ];do


		if [ $mode == "somatic" ];then
	        	mkdir tmp_tumor;
		        mkdir tmp_normal;
		        tmp_folder="$PWD/tmp/"
		        tmp_folder1="$PWD/tmp_tumor/"
		        tmp_folder2="$PWD/tmp_normal/"
		        Rscript $library/SV_local_CN_segment_wgs_raw_somatic.R $SV $i $bin $tmp_folder $bicseq $lambda_ini $bam $normal_bin $tmp_folder1 $tmp_folder2 $normal_bam $chr_prefix
			rm -rf tmp_tumor;
			rm -rf tmp_normal;
		else
			tmp_folder="$PWD/tmp/"
			Rscript $library/SV_local_CN_segment_wgs_raw.R $SV $i $bin $tmp_folder $bicseq $lambda_ini $bam $chr_prefix;
		fi
	
		abs_line=`cat copy_numbers | wc -l`
		if [ $abs_line -gt 50000 ];then
			lambda_ini=$(($lambda_ini*2));
			echo $lambda_ini >> "lambda.log";
		fi
		rm -rf tmp;
		
	done


	Rscript $library/ABSOLUTE.R "$cancer_type" $sample $min_ploidy $max_ploidy $ABSOLUTE_STABLE
	if [[ -s ABSOLUTE_output/output/reviewed/$sample.test.ABSOLUTE.table.txt ]];then
		echo "pass"
	else
		ABSOLUTE_STABLE="T";
	        Rscript $library/ABSOLUTE.R "$cancer_type" $sample $min_ploidy $max_ploidy $ABSOLUTE_STABLE
	        if [[ -s ABSOLUTE_output/output/reviewed/$sample.test.ABSOLUTE.table.txt ]];then
			echo "pass"
		else
			exit 1
		fi
	fi
	
	
	Rscript $library/integer_programming.R $sample "copy_numbers_ABSOLUTE_input.negative_marker" $SV $read_length> IP.log
	test1=`cat SVs.CN_opt.filtered | wc -l`
	test2=`cat SVs | wc -l`
	
	cp SVs.CN_opt.filtered ../SVs
	cd ../
	if [ $test1 -eq $test2 ];then
		break;
	fi
done


cd iter$iter

	Rscript $library/remapping/unbalanced_nodes.R
	perl $library/remapping/edge_filling_job_optimized.pl $fasta_prefix.fa $NPE_dir $haplotype_coverage $library $search_length $CIGAR $read_length $chr_prefix> edge_probability
	perl $library/remapping/edge_filling_job_blat_recal.pl $fasta_prefix
	perl $library/remapping/edge_filling_job_blat_recal_sam_To_bed_re.pl $fasta_prefix.fa $CIGAR $read_length $search_length > edge_filling_sorted_paired_primary.sam.info.recal
	perl $library/remapping/edge_filling_job_optimized_for_recal_info.pl $library $haplotype_coverage $search_length > edge_probability_recal
	sort -nrk 5,5 edge_probability_recal  > edge_probability_recal_sorted
	cat  edge_probability_recal_sorted | awk '{if($5>0.01) print $0}' > edge_probability_recal_sorted.filtered
	cp SVs SVs.added
	perl $library/remapping/module_edge_to_SVs.pl $search_length >> SVs.added
	cp SVs.added ../SVs
cd ../


iter_start=$iter+1;

for (( iter=$iter_start; iter<=100; iter++ ))
do
        mkdir iter$iter;
        cd iter$iter;
        cp  ../SVs ./
        mkdir tmp;

        if [ $mode == "somatic" ];then
                mkdir tmp_tumor;
                mkdir tmp_normal;
                tmp_folder="$PWD/tmp/"
                tmp_folder1="$PWD/tmp_tumor/"
                tmp_folder2="$PWD/tmp_normal/"
                Rscript $library/SV_local_CN_segment_wgs_raw_somatic.R $SV $i $bin $tmp_folder $bicseq $lambda_fi $bam $normal_bin $tmp_folder1 $tmp_folder2 $normal_bam $chr_prefix
		rm -rf tmp_tumor;
		rm -rf tmp_normal;

        else
                tmp_folder="$PWD/tmp/"
                Rscript $library/SV_local_CN_segment_wgs_raw.R $SV $i $bin $tmp_folder $bicseq $lambda_fi $bam $chr_prefix;
        fi
	rm -rf tmp;

        Rscript $library/ABSOLUTE.R $cancer_type $sample $min_ploidy $max_ploidy $ABSOLUTE_STABLE
        if [[ -s ABSOLUTE_output/output/reviewed/sample1.test.ABSOLUTE.table.txt ]];then
                echo "pass"
        else
                ABSOLUTE_STABLE="T";
                Rscript $library/ABSOLUTE.R $cancer_type $sample $min_ploidy $max_ploidy $ABSOLUTE_STABLE
                if [[ -s ABSOLUTE_output/output/reviewed/sample1.test.ABSOLUTE.table.txt ]];then
                        echo "pass"
                else
                        exit 1
                fi
        fi

        Rscript $library/integer_programming.R $sample "copy_numbers_ABSOLUTE_input.negative_marker" $SV $read_length > IP.log
        test1=`cat SVs.CN_opt.filtered | wc -l`
        test2=`cat SVs | wc -l`

        cp SVs.CN_opt.filtered ../SVs
        cd ../
        if [ $test1 -eq $test2 ];then
                break;
        fi
done
