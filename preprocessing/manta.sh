#!/bin/bash
script1=configManta.py
tumor_bam=$1
reference=$2
PE_thres=4

mkdir -p log
$script1 --tumorBam $tumor_bam --referenceFasta $reference --runDir manta &> log/manta.log
cd manta
./runWorkflow.py -m local -j 8 &>> ../log/manta.log
cd ../


zcat ./manta/results/variants/tumorSV.vcf.gz | awk -F "\t" '{

	n=split($8,f,";");
	split(f[2],DUPDEL,"=");
	split(f[1],TRA, "=");
	rn=split($10,read, ":");
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

}' | awk '{ if($7 > '$PE_thres' ){
        print $0"\t0\t0\t0\t0\t0\t0"}}' | awk '{if(($2~/^[1-2]?[0-9]$/ || $2=="X") && ($4~/^[1-2]?[0-9]$/ || $4=="X")) print $0}' > manta.format
