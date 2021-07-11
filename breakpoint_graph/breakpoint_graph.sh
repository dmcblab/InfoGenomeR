#!/bin/bash
lambda_ini=1
lambda_fi=4
search_length=50000;
CIGAR="100M";
read_length=100;
min_ploidy=1.5
max_ploidy=5
sample=sample1
cancer_type=null
npe_dir=null
fasta_prefix="hg19"
out_dir=InfoGenomeR_job
chr_prefix=0  # The reference genome has "chr" prefix or not

bam=null
normal_bam=null
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -m|--mode)
      mode="$2"
      shift
      shift
      ;;
    -i|--lambda_ini)
      lambda_ini="$2"
      shift
      shift
      ;;
    -f|--lambda_fi)
      lambda_fi="$2"
      shift
      shift
      ;;
    -t|--cancer_type)
      cancer_type="$2"
      shift
      shift
      ;;
    -d|--npe_dir)
      npe_dir=`readlink -f "$2"`
      shift
      shift
      ;;
    -g|--ref_genome)
      fasta_prefix=`readlink -f "$2"`
      shift
      shift
      ;;
    -c|--cn_norm_germ)
      cn_norm_germ=`readlink -f "$2"`
      shift
      shift
      ;;
    -s|--germ_LocSeq_result)
      germ_LocSeq_result=`readlink -f "$2"`
      shift
      shift
      ;;
    -o|--out_dir)
      out_dir="$2"
      shift
      shift
      ;;
    -h|--help)
      echo -e "Usage: ./breakpoint_graph.sh <SVs> <cn_norm> [options]\n"
      echo -e "Options:"
      echo -e  "\t-m, --mode (required)\n \t\t Select the mode (germline, total, somatic)
\t-i, --lambda_ini (required) \n \t\t Initial lambda for the first round iterations (default: 1)
\t-f, --lambda_fi (required) \n \t\t Final lambda for the second round iterations (default: 4)
\t-t, --cancer_type (optional) \n \t\t Cancer type used for ABSOLUTE estimation (BRCA, GBM, OV, ...). If unknown, write null.
\t-d, --npe_dir (optional) \n \t\t Directory that contains NPE.fq1 and NPE.fq2 (non-properly paired reads). If it is not assigned, InfoGenomeR runs without NP reads mapping.
\t-g, --ref_genome (required for NP reads mapping) \n \t\t Fasta prefix (hg19 or hg38). Enter the prefix without .2bit and .fa extension.
\t-c, cn_norm_germ (required for somatic mode) \n \t\t Directory that contains copy number bins from a control genome.
\t-s, --germ_LocSeq_result (required for somatic mode) \n \t\t Local segmentation results from a control genome.
\t-o, --out_dir \n \t\t If it already exists, results are overlaid (default: InfoGenomeR_job)
\t-h, --help\n
"
      
      exit 1
      ;;
    *)  
     POSITIONAL+=("$1") # save it in an array for later
     shift
     ;;
  esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters

echo "InfoGenomeR_lib = ${InfoGenomeR_lib}"
if [[ ! -s ${InfoGenomeR_lib}/breakpoint_graph ]];then
	echo "set the InfoGenomeR_lib path correctly"
	exit 1
fi
echo "BICseq2_path = ${BICseq2_path}"
if [[ ${BICseq2_path} == "" ]];then
        echo "set the InfoGenomeR_lib path"
        exit 1
fi
echo "mode = ${mode}"
if [[ ${mode} == "" ]];then
	echo "set the mode; total, germline, or somatic"
	exit 1
fi
echo "lambda_ini  = ${lambda_ini}"
echo "lambda_fi   = ${lambda_fi}"
echo "cancer_type = ${cancer_type}"
if [[ -s $npe_dir\/NPE.fq1 ]] && [[ -s $npe_dir\/NPE.fq2 ]];then
	echo "NP reads mapping = True"
else
	npe_dir="null"
	echo "NPE_mapping = False"
fi
echo "NPE_dir = ${npe_dir}"
echo "fasta_prefix = ${fasta_prefix}"
if [[ -s $npe_dir ]];then
	if [[ ! -s ${fasta_prefix}.fa ]] || [[ ! -s ${fasta_prefix}.2bit ]];then
		echo "write the reference genome prefix and check .2bit exists"
		exit 1
	fi
fi
SV=`readlink -f "$1"`
echo -e "SVs = $SV"
if [[ $1 == "" ]];then
	echo "An SV file is required"
	exit 1
fi
cn_norm=`readlink -f "$2"`
echo "cn_norm = $cn_norm"
if [[ $2 == "" ]];then
        echo "cn_norm directory is required"
        exit 1
fi
echo "cn_norm_germ = ${cn_norm_germ}"
echo "germ_LocSeq_result = ${germ_LocSeq_result}"
if [[ $mode == "somatic" ]];then
	if [[ $cn_norm_germ == "" ]];then
		echo "cn_norm_germ directory is required";
		exit 1
	fi
        if [[ $germ_LocSeq_result == "" ]];then
                echo "germ_LocSeq_result is required";
                exit 1
        fi
fi

bin="$cn_norm/"
normal_bin="$cn_norm_germ/"
if [[ -s $out_dir ]];then
	echo "out_dir exists... overlapping"
fi
mkdir -p $out_dir
cd $out_dir
cp $SV ./SVs
SV=SVs
i=1;
iter=1;
library=$InfoGenomeR_lib\/breakpoint_graph
bicseq=$BICseq2_path\/NBICseq-seg.pl
if [ $mode == "germline" ];then
	lambda_ini=1;
	Rscript $library/SV_break_adjustment.R
	$library/SV_truncated_wgs.sh SVs $bin
	Rscript $library/SV_break_adjust.R
	mkdir -p tmp;
	tmp_folder="$PWD/tmp/"
        Rscript $library/SV_local_CN_segment_wgs_raw.R $SV $i $bin $tmp_folder $bicseq $lambda_ini $bam $chr_prefix;
	rm -rf tmp;
	exit 1
fi



if [ $mode == "somatic" ];then
	mkdir -p cn_norm
	mkdir -p cn_norm_germ
	echo "filtering germline bins"
	for chr in {1..23}; do
		Rscript $library/germ_cnv_filter_v2.R $germ_LocSeq_result $normal_bin $bin $chr
		echo "$chr.bin"
	done
	
	bin="$PWD/cn_norm/"
	normal_bin="$PWD/cn_norm_germ/"
fi


Rscript $library/SV_break_adjustment.R
$library/SV_truncated_wgs.sh SVs $bin


ABSOLUTE_STABLE="F";

ROUND=1
for iter in {1..100}
do
	mkdir -p iter$iter;
	cd iter$iter;
	cp  ../SVs ./
	mkdir -p tmp;

	if [[ $ROUND -eq 1 ]];then
		lambda=$lambda_ini
	else
		lambda=$lambda_fi
	fi

	echo "Local segmentation..."
	abs_line=50000;
        while [ $abs_line -ge 50000 ];do


		if [ $mode == "somatic" ];then
	        	mkdir -p tmp_tumor;
		        mkdir -p tmp_normal;
		        tmp_folder="$PWD/tmp/"
		        tmp_folder1="$PWD/tmp_tumor/"
		        tmp_folder2="$PWD/tmp_normal/"
		        Rscript $library/SV_local_CN_segment_wgs_raw_somatic.R $SV $i $bin $tmp_folder $bicseq $lambda $bam $normal_bin $tmp_folder1 $tmp_folder2 $normal_bam $chr_prefix
			rm -rf tmp_tumor;
			rm -rf tmp_normal;
		else
			tmp_folder="$PWD/tmp/"
			Rscript $library/SV_local_CN_segment_wgs_raw.R $SV $i $bin $tmp_folder $bicseq $lambda $bam $chr_prefix;
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
	if [[ $test1 -eq $test2 ]] && [[ $ROUND -eq 1 ]];then
		if [[ -s $npe_dir ]];then
			cd iter$iter
			total_coverage_h=`cat $bin/*.norm.bin  | awk 'BEGIN{sum=0;n=0;}{sum=sum+($2-$1)*$3/100; n=n+1}END{print sum/n}'`
			ploidy_h=`cat ABSOLUTE_output/output/reviewed/*.test.ABSOLUTE.table.txt | awk -F "\t" '{print $5}' | tail -n 1`
			purity_h=`cat ABSOLUTE_output/output/reviewed/*.test.ABSOLUTE.table.txt | awk -F "\t" '{print $4}' | tail -n 1`
			bc_h=`echo -e "$total_coverage_h\t$ploidy_h\t$purity_h" | awk '{b=$1/($2*$3+2*(1-$3)); print (b-int(b)<0.499)?int(b):int(b)+1}'`
			base_coverage=$bc_h
			Rscript $library/remapping/unbalanced_nodes.R
			perl $library/remapping/edge_filling_job_optimized.pl $fasta_prefix.fa $npe_dir $base_coverage $library $search_length $CIGAR $read_length $chr_prefix> edge_probability
			perl $library/remapping/edge_filling_job_blat_recal.pl $fasta_prefix
			perl $library/remapping/edge_filling_job_blat_recal_sam_To_bed_re.pl $fasta_prefix.fa $CIGAR $read_length $search_length > edge_filling_sorted_paired_primary.sam.info.recal
			perl $library/remapping/edge_filling_job_optimized_for_recal_info.pl $library $base_coverage $search_length > edge_probability_recal
			sort -nrk 5,5 edge_probability_recal  > edge_probability_recal_sorted
			cat  edge_probability_recal_sorted | awk '{if($5>0.01) print $0}' > edge_probability_recal_sorted.filtered
			cp SVs SVs.added
			perl $library/remapping/module_edge_to_SVs.pl $search_length >> SVs.added
			cp SVs.added ../SVs
			cd ../
		fi
		if [[ ! -s iter$iter\/SVs.added ]];then
                        cd iter$iter
			Rscript $library/TO_bpgraph_hidden_node_for_connectivity.R
			Rscript $library/lite_SV_edges.R ../iter1/SVs $search_length F
			cp SVs.added ../SVs
			cd ../

		fi
		ROUND=2
	elif [[ $test1 -eq $test2 ]];then
                echo "Breakpoint graph construction is finished"
		break;
	fi
done
