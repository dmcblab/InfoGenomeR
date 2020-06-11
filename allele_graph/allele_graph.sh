#!/bin/bash
library="/home/qlalf1457/InfoGenomer/allele_graph"
snp_count_path="/home/qlalf1457/beagle_1000_Genomes_phase3"

mode=$1;
germ_LocSeq_result=$2;
hom=$3;
het=$4;
ref=$5;

iter=`ls -l | grep -E 'iter[1-9]?[0-9]?[0-9]$' | awk 'BEGIN{max=0}{split($9,f,"iter"); if(max<f[2]) max=f[2];}END{print max}'`

rm exclude
echo -n "" > exclude
if [[ $mode == "somatic" ]];then
        cat $germ_LocSeq_result | awk '{if($6 < -0.5 || $6 > 0.3 ) print $2"\t"$3"\t"$4}' > exclude
fi

final_iter=$(($iter -1 ));

for i in `seq 1 $final_iter`;do
        if [[ -s iter$i\/exclude ]];then
                cat iter$i\/exclude >> exclude
        fi
done

cd iter$iter
$library/HR_tag.sh $library $ref;

cp ../exclude ./
cp $hom ./hom_snps.format
cp $het ./het_snps.format
if [[ -s exclude ]];then
	Rscript $library/snps_remove.R
	$library/ACN_estimation.sh $library T $snp_count_path
else
        $library/ACN_estimation.sh $library F $snp_count_path
fi

Rscript $library/SNP_phasing_decision_boundary_ver_with_homozygous.R $snp_count_path
Rscript $library/genotype_format_for_popul_phasing.R  $snp_count_path

