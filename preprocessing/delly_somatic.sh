#!/bin/bash
lib=`dirname $(readlink -f ${BASH_SOURCE[0]})`

cleanup() {
	pkill -P $$
	kill 0
}

for sig in INT QUIT HUP TERM; do
	trap "
	    cleanup
            trap - $sig EXIT
            kill -s $sig "'"$$"' "$sig"
done

# trap cleanup EXIT

tumor_bam=$1
normal_bam=$2
ref=$3
PE_thres=$4
MAPQ_thres=$5

idx=0;

declare -a pids;
mkdir -p log
for type in DEL INV DUP TRA;do
	delly call -t $type -g $ref -o somatic_pre_$type.bcf $tumor_bam $normal_bam &>log/somatic.$type.log &
	pids[$idx]=$!;
	idx=$(($idx+1));
done

for z in `seq 0 3`;do
        wait ${pids[$z]};
done


bcftools view somatic_pre_DEL.bcf -h | tail -n1 | awk '{print $10"\ttumor"; print $11"\tcontrol"}' > samples.tsv

delly filter -t DEL -f somatic -o somatic_DEL.bcf -s samples.tsv somatic_pre_DEL.bcf -a 0 -m 100
delly filter -t INV -f somatic -o somatic_INV.bcf -s samples.tsv somatic_pre_INV.bcf -a 0 -m 100
delly filter -t DUP -f somatic -o somatic_DUP.bcf -s samples.tsv somatic_pre_DUP.bcf -a 0 -m 100
delly filter -t TRA -f somatic -o somatic_TRA.bcf -s samples.tsv somatic_pre_TRA.bcf -a 0 -m 100

for type in DEL INV DUP TRA;do
	bcftools view somatic_$type.bcf > somatic_$type.vcf
done

echo -n "" > delly.format
for type in DEL INV DUP TRA;
do
                tumorID=somatic_$type.vcf

                cat $tumorID | grep -v '#' | awk -F "\t" '{
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

Rscript $lib/SV_overnoise.R delly.format


