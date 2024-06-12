#!/bin/bash
unset R_HOME
mode=$1
root_dir=`readlink -f $2`
preprocess_dir=`readlink -f $3`
output_dir=`readlink -f $4`

ref=`readlink -f $5`
ref_dir=`dirname $ref`
ref_prefix=`basename $ref | awk -F "." '{print $1}'`

lambda_ini=$6
lambda_fi=$7
min_ploidy=$8
max_ploidy=$9
cancer_type="${10}"

if [[ $ref =~ "hg19" ]];then
	export Ref_version=GRCh37
fi

bicseq_script=`which NBICseq-norm.pl`
bicseq_dir=`dirname $bicseq_script`

LIB=`readlink -f ${BASH_SOURCE[0]} | awk '{n=split($1,f,"/"); for(i=1;i<=n-3;i++){printf "%s/", f[i]}}'`
export BICseq2_path=$bicseq_dir
export InfoGenomeR_lib=$LIB
export PATH=$InfoGenomeR_lib/breakpoint_graph:$InfoGenomeR_lib/allele_graph:$InfoGenomeR_lib/haplotype_graph:$InfoGenomeR_lib/Eulerian:$PATH


cd $root_dir
mkdir -p log

if [[ $mode == "somatic" ]];then
	breakpoint_graph -m germline ${preprocess_dir}/delly/germline.filtered.format  -o germline_job ${preprocess_dir}/bicseq/cn_norm_germ &> log/germline_job.log
	breakpoint_graph -m somatic ${preprocess_dir}\/SVs ${preprocess_dir}/bicseq/cn_norm -t "$cancer_type" -i $lambda_ini -f $lambda_fi -n $min_ploidy -x $max_ploidy -g $ref_dir\/$ref_prefix  -c ${preprocess_dir}/bicseq/cn_norm_germ -s germline_job/copy_numbers -o InfoGenomeR_job &> log/InfoGenomeR_job.log
else
	breakpoint_graph -m total ${preprocess_dir}\/SVs ${preprocess_dir}/bicseq/cn_norm -t "$cancer_type" -i $lambda_ini -f $lambda_fi -n $min_ploidy -x $max_ploidy -g $ref_dir\/$ref_prefix -o InfoGenomeR_job &> log/InfoGenomeR_job.log
fi


