#!/bin/bash
reference=hg19.fa
script1=./manta-1.1.0.centos5_x86_64/bin/configManta.py
script2=./runWorkflow.py
normal_bam=normal.bam
tumor_bam=tumor.bam

$script1 --normalBam $normal_bam --tumorBam $tumor_bam --referenceFasta $reference --runDir manta_somatic
$script2 -m local -j 1 


zcat ./manta_somatic/results/variants/somaticSV.vcf.gz | awk -F "\t" '{
	n=split($8,f,";");
	split(f[2],DUPDEL,"=");
	split(f[1],TRA, "=");
	rn=split($11,read, ":");
	if(rn==1){
		split(read[1],pr,",");
		PR=pr[2];
		SR=0;
	}else{
		split(read[1],pr,",");
		PR=pr[2];
		split(read[2],sr,",");
		SR=sr[2];
	}
	if(DUPDEL[2]=="DUP"){
		split(f[1], pos2, "=");
		print "<DUP>\t"$1"\t"$2"\t"$1"\t"pos2[2]"\t5to3\t"PR"\t"SR
	}else if(DUPDEL[2]=="DEL"){
		split(f[1], pos2, "=");
		print "<DEL>\t"$1"\t"$2"\t"$1"\t"pos2[2]"\t3to5\t"PR"\t"SR
	}else if(DUPDEL[2]=="INV"){
		split(f[1], pos2, "=");
		split($8, INV, ";INV");
		split(INV[2], inv_ori, ";");
		print "<INV>\t"$1"\t"$2"\t"$1"\t"pos2[2]"\t"inv_ori[1]"to"inv_ori[1]"\t"PR"\t"SR;
	}else if(TRA[2]=="BND"){
		n1=split($5, tra_ori, "[");
		if(n1>2){
			if(length(tra_ori[1])!=0)
				ori="3to5"
			else if(length(tra_ori[3])!=0)
				ori="5to5"
		}else{
		n2=split($5, tra_ori, "]");
			if(length(tra_ori[1])!=0)
				ori="3to3"
			else if(length(tra_ori[3])!=0)
				ori="5to3"
		}
		split(tra_ori[2], chr2_pos2, ":");
		print "<TRA>\t"$1"\t"$2"\t"chr2_pos2[1]"\t"chr2_pos2[2]"\t"ori"\t"PR"\t"SR;


	}

}' | awk '{
	print $0"\t0\t0\t0\t0\t0\t0"}'  > manta.format

