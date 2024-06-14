#!/bin/bash
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

reference=$1
normal_bam=$2
tumor_bam=$3

declare -a pids;

mkdir -p vcf_somatic
for l in `seq 1 22`;do
		samtools mpileup -t DP,AD,ADF,ADR,SP -q 10 -d 5000 -I -uf $reference $normal_bam $tumor_bam -r $l | bcftools call -c -v > vcf_somatic/$l.vcf &
		pids[$l]=$!;
done
l='X';
                samtools mpileup -t DP,AD,ADF,ADR,SP -q 10 -d 5000 -I -uf $reference $normal_bam $tumor_bam -r $l | bcftools call -c -v > vcf_somatic/$l.vcf &
		pids[23]=$!;

for l in `seq 1 23`;do
	wait ${pids[$l]}
done


for i in `seq 1 23`; do
        if [ $i -eq 23 ];then
                i="X";
        fi
        cat vcf_somatic/$i.vcf | grep -v '#' | awk '{

                if(length($4)==1 && length($5) == 1 && $6>50){
                        split($10,f,":");
                        split($11,g,":");
                        if(f[1] == "0/1"){
                                split($10,f1,":");
                                split(f1[7], normal, ",");
                               split($11,f2,":");
                                split(f2[7], tumor, ",");

                                print "'$i'""\t"$2"\t"$3"\t"$4"\t"$5"\t"normal[1]"\t"normal[2]"\t"tumor[1]"\t"tumor[2]
                        }
                }
        }' > vcf_somatic/$i.vcf.format.het

        cat vcf_somatic/$i.vcf | grep -v '#' | awk '{

                if(length($4)==1 && length($5) == 1 && $6>50){
                        split($11,f,":");
                        if(f[1] == "1/1"){
                                split($10,f1,":");
                                split(f1[7], normal, ",");
                              split($11,f2,":");
                                split(f2[7], tumor, ",");

                                print "'$i'""\t"$2"\t"$3"\t"$4"\t"$5"\t"normal[1]"\t"normal[2]"\t"tumor[1]"\t"tumor[2]
                        }
                }
        }' > vcf_somatic/$i.vcf.format.hom

done



cat vcf_somatic/*.vcf.format.het > het_snps.format
cat vcf_somatic/*.vcf.format.hom > hom_snps.format

