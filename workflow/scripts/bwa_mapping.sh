#!/bin/bash
declare -a pids;

fastq_dir=`readlink -f $1`
REF=`readlink -f $2`
out_dir=`readlink -f $3`

tumor_fq1=$fastq_dir\/tumor1.fq.gz
tumor_fq2=$fastq_dir\/tumor2.fq.gz

normal_fq1=$fastq_dir\/normal1.fq.gz
normal_fq2=$fastq_dir\/normal2.fq.gz


mkdir -p $out_dir



bwa mem -t 40 $REF $tumor_fq1 $tumor_fq2 | samtools view -bS > $out_dir\/tumor.bam
samtools sort $out_dir\/tumor.bam -@ 40 > $out_dir\/tumor_sorted.bam
samtools index $out_dir\/tumor_sorted.bam


bwa mem -t 40 $REF $normal_fq1 $normal_fq2 | samtools view -bS > $out_dir\/normal.bam
samtools sort $out_dir\/normal.bam -@ 40 > $out_dir\/normal_sorted.bam
samtools index $out_dir\/normal_sorted.bam 


rm $out_dir\/tumor.bam
rm $out_dir\/normal.bam
