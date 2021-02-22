#!/bin/bash

#echo "Type SRR ID"
#read tumorIDini
ref=hg19.fa
tsv=samples.tsv
tumor_bam=tumor.bam
normal_bam=normal.bam

delly_v0.7.6 call -t DEL -g $ref -o somatic_pre_DEL.bcf $tumor_bam $normal_bam
delly_v0.7.6 filter -t DEL -f somatic -o somatic_DEL.bcf -s samples.tsv somatic_pre_DEL.bcf -a 0 -m 100
delly_v0.7.6 call -t INV -g $ref -o somatic_pre_INV.bcf $tumor_bam $normal_bam
delly_v0.7.6 filter -t INV -f somatic -o somatic_INV.bcf -s samples.tsv somatic_pre_INV.bcf -a 0 -m 100
delly_v0.7.6 call -t DUP -g $ref -o somatic_pre_DUP.bcf $tumor_bam $normal_bam
delly_v0.7.6 filter -t DUP -f somatic -o somatic_DUP.bcf -s samples.tsv somatic_pre_DUP.bcf -a 0 -m 100
delly_v0.7.6 call -t TRA -g $ref -o somatic_pre_TRA.bcf $tumor_bam $normal_bam 
delly_v0.7.6 filter -t TRA -f somatic -o somatic_TRA.bcf -s samples.tsv somatic_pre_TRA.bcf -a 0 -m 100

PE_thres=$1
MAPQ_thres=$2
#lib=$2

for i in {1..4}
do
        if [ $i == 1 ]
        then
                tumorID=somatic_DEL.vcf
        elif [ $i == 2 ]
        then
                tumorID=somatic_INV.vcf
        elif [ $i == 3 ]
        then
                tumorID=somatic_DUP.vcf
        else
                tumorID=somatic_TRA.vcf
        fi

                cat $tumorID | awk -F "\t" '{
                                                n=split($8,f,"SR=");
                                                if(n>1){
                                                        split(f[2],SR, ";");
                                                        SR_number=SR[1];
                                                }else{
                                                        SR_number=0;
                                                }
                                                split($8, f, "SVTYPE=");
                                                split(f[2], SVTYPE, ";");

                                               split($8, f, "CHR2=");
                                                split(f[2], chr2, ";");

                                               split($8, f, ";END=");
                                                split(f[2], pos2, ";");

                                               split($8, f, ";PE=");
                                                split(f[2], PE, ";");

                                               split($8, f, "MAPQ=");
                                                split(f[2], MAPQ, ";");

                                              split($8, f, "CT=");
                                                split(f[2], ori, ";");

                                                if($7=="PASS" && (PE[1] > '$PE_thres' && SR_number > -1)&& MAPQ[1] > '$MAPQ_thres'){
                                                        print "<"SVTYPE[1]">""\t"$1"\t"$2"\t"chr2[1]"\t"pos2[1]"\t"ori[1]"\t"PE[1]"\t"SR_number"\t"MAPQ[1]"\t"PE[1]+SR_number"\t0\t0\t1\t1\t"
                                                }
                                        }' | awk '{if(($2~/^[1-2]?[0-9]$/ || $2=="X") && ($4~/^[1-2]?[0-9]$/ || $4=="X")) print $0}'  >> delly.format

done

#Rscript $lib/SV_overnoise.R delly.format


