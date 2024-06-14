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


mode=$1
bam=$2
ref=$3
exclude_tsv=$4
PE_thres=$5
MAPQ_thres=$6

bcf_prefix="";
if [[ $mode == "germline" ]];then
	bcf_prefix="germline"
else
	bcf_prefix="tumor"
fi

idx=0;
mkdir -p log;

for type in DEL DUP INV TRA;do
	delly call -t $type -g $ref -o $bcf_prefix.$type.bcf $bam  -x $exclude_tsv -q 10 -s 15 -n &>log/$bcf_prefix.$type.log &
	pids[$idx]=$!;
	idx=$(($idx+1));
done


for pid in ${pids[*]}; do
        wait $pid;
done

for type in DEL DUP INV TRA;do
	bcftools view $bcf_prefix.$type.bcf > $bcf_prefix.$type.vcf
done

echo -n "" > $bcf_prefix.filtered.format
for type in DEL DUP INV TRA;do
                file=$bcf_prefix.$type.vcf

                cat $file | awk -F "\t" '{
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
                                        }' | awk '{if(($2~/^[1-2]?[0-9]$/ || $2=="X") && ($4~/^[1-2]?[0-9]$/ || $4=="X")) print $0}'  >> $bcf_prefix.filtered.format

done

if [[ mode != "germline" ]];then
	ln -s $bcf_prefix.filtered.format delly.format
fi
