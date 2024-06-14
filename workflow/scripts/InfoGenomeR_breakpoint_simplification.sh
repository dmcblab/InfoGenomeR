#!/bin/bash
unset R_HOME
mode=$1
root_dir=`readlink -f $2`
InfoGenomeR_dir=`readlink -f $3`
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
iter=`ls ${InfoGenomeR_dir} -l | grep -E 'iter[1-9]?[0-9]?[0-9]$' | awk 'BEGIN{max=0}{split($9,f,"iter"); if(max<f[2]) max=f[2];}END{print max}'`
mkdir -p InfoGenomeR_simplification_job


cp -r ${InfoGenomeR_dir}/iter$iter InfoGenomeR_simplification_job
cp ${InfoGenomeR_dir}/SVs InfoGenomeR_simplification_job
ln -s ${InfoGenomeR_dir}/cn_norm InfoGenomeR_simplification_job/cn_norm


if [[ $mode == "somatic" ]];then
	ln -s ${InfoGenomeR_dir}/cn_norm_germ InfoGenomeR_simplification_job/cn_norm_germ
	breakpoint_graph -m simplification -o InfoGenomeR_simplification_job  -t "$cancer_type" -i $lambda_ini -f $lambda_fi -n $min_ploidy -x $max_ploidy InfoGenomeR_simplification_job/SVs InfoGenomeR_simplification_job/cn_norm -g $ref_dir\/$ref_prefix -c InfoGenomeR_simplification_job/cn_norm_germ -s germline_job/copy_numbers &>  log/simplification.log
else
       breakpoint_graph -m simplification -o InfoGenomeR_simplification_job  -t "$cancer_type" -i $lambda_ini -f $lambda_fi -n $min_ploidy -x $max_ploidy InfoGenomeR_simplification_job/SVs InfoGenomeR_simplification_job/cn_norm -g $ref_dir\/$ref_prefix  &>  log/simplification.log
fi


