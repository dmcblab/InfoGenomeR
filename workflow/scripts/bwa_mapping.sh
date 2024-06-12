#!/bin/bash
unset R_HOME
declare -a pids;

mode=$1
fastq_dir=`readlink -f $2`
REF=`readlink -f $3`
out_dir=`readlink -f $4`


if [[ $mode != "somatic" ]] && [[ $mode != "total" ]];then
	echo "mode should be somatic or total"
	exit 1
fi


tumor_fq1=$fastq_dir\/tumor1.fq.gz
tumor_fq2=$fastq_dir\/tumor2.fq.gz

if [[ $mode == "somatic" ]];then
	normal_fq1=$fastq_dir\/normal1.fq.gz
	normal_fq2=$fastq_dir\/normal2.fq.gz
fi

mkdir -p $out_dir

bwa mem -t 40 $REF $tumor_fq1 $tumor_fq2 | samtools view -bS > $out_dir\/tumor.bam
samtools sort $out_dir\/tumor.bam -@ 40 > $out_dir\/tumor_sorted.bam
samtools index $out_dir\/tumor_sorted.bam
rm $out_dir\/tumor.bam

if [[ $mode == "total" ]];then
  echo "tumor read mapping finished"
  exit 0
fi

bwa mem -t 40 $REF $normal_fq1 $normal_fq2 | samtools view -bS > $out_dir\/normal.bam
samtools sort $out_dir\/normal.bam -@ 40 > $out_dir\/normal_sorted.bam
samtools index $out_dir\/normal_sorted.bam 

rm $out_dir\/normal.bam

echo "normal read mapping finished"
