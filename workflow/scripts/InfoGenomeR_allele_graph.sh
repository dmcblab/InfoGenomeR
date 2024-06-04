#!/bin/bash

root_dir=`readlink -f $1`
InfoGenomeR_dir=`readlink -f $2`
output_dir=`readlink -f $3`

ref=`readlink -f $4`
ref_dir=`dirname $ref`
ref_prefix=`basename $ref | awk -F "." '{print $1}'`

export Haplotype_path=`readlink -f $5`


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


mkdir -p InfoGenomeR_allele_graph_job

cp -r ${InfoGenomeR_dir}/iter$iter InfoGenomeR_allele_graph_job

allele_graph -m somatic -g $ref_dir\/$ref_prefix -o InfoGenomeR_allele_graph_job -t 23 -s germline_job/copy_numbers preprocess/snp/hom_snps.format preprocess/snp/het_snps.format &> log/allele_graph.log

