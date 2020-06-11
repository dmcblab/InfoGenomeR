#!/bin/bash

library="/home/qlalf1457/InfoGenomer/haplotype_graph/phasing"
snp_count_path="/home/qlalf1457/beagle_1000_Genomes_phase3"

iter=`ls -l | grep -E 'iter[1-9]?[0-9]?[0-9]$' | awk 'BEGIN{max=0}{split($9,f,"iter"); if(max<f[2]) max=f[2];}END{print max}'`

cd iter$iter


for i in `seq 1 22`;do
	$library/phasing $snp_count_path/chr$i.chr$i.bgl.dag.align genotype.format.$i  >  genotype.format.$i.phased && paste $snp_count_path/chr$i.marker   genotype.format.$i.phased  | cut -f1,2,4,6,7 > diploid_haplotype.$i
done


i="X"
        $library/phasing $snp_count_path/chr$i.chr$i.bgl.dag.align genotype.format.$i  >  genotype.format.$i.phased && paste $snp_count_path/chr$i.marker   genotype.format.$i.phased  | cut -f1,2,4,6,7 > diploid_haplotype.$i



