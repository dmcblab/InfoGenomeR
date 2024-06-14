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
tumor_bam=$2

declare -a pids;

mkdir -p vcf

for l in `seq 1 22`;do
		samtools mpileup -t DP,AD,ADF,ADR,SP -q 20 -d 5000 -I -uf $reference $tumor_bam -r $l | bcftools call -c -v > vcf/$l.vcf &
		pids[$l]=$!;
done
l='X';
                samtools mpileup -t DP,AD,ADF,ADR,SP -q 20 -d 5000 -I -uf $reference $tumor_bam -r $l | bcftools call -c -v > vcf/$l.vcf &
		pids[23]=$!;

for l in `seq 1 23`;do
	wait ${pids[$l]}
done


for i in `seq 1 23`; do
        if [ $i -eq 23 ];then
                i="X";
        fi
        cat vcf/$i.vcf | grep -v '#' | awk '{

                if(length($4)==1 && length($5) == 1 && $6>50){
                        split($10,f,":");
                        if(f[1] == "0/1"){
                                split($10,f1,":");
                                split(f1[7], tumor, ",");

                                print "'$i'""\t"$2"\t"$3"\t"$4"\t"$5"\t""NA""\t""NA""\t"tumor[1]"\t"tumor[2]
                        }
                }
        }' > vcf/$i.vcf.format.het

        cat vcf/$i.vcf | grep -v '#' | awk '{

                if(length($4)==1 && length($5) == 1 && $6>50){
                        split($10,f,":");
                        if(f[1] == "1/1"){
                                split($10,f1,":");
                                split(f1[7], tumor, ",");

                                print "'$i'""\t"$2"\t"$3"\t"$4"\t"$5"\t""NA""\t""NA""\t"tumor[1]"\t"tumor[2]
                        }
                }
        }' > vcf/$i.vcf.format.hom


done

cat vcf/*.vcf.format.het > het_snps.format
cat vcf/*.vcf.format.hom > hom_snps.format

