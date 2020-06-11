#!/bin/bash
library="/home/qlalf1457/InfoGenomer/haplotype_graph"
snp_count_path="/home/qlalf1457/beagle_1000_Genomes_phase3"

iter=`ls -l | grep -E 'iter[1-9]?[0-9]?[0-9]$' | awk 'BEGIN{max=0}{split($9,f,"iter"); if(max<f[2]) max=f[2];}END{print max}'`

./phasing/phasing.sh

cd iter$iter

Rscript $library/ACN_order_by_haplotype.R
Rscript $library/AS_SV_haplotype_phased.R
#Rscript $library/ACN_GRAPH_phased_ver_haplotype_phased.R
#Rscript $library/ACN_GRAPH_phased_ver_haplotype_phased_hilight.R F
Rscript $library/TO_bpgraph_hidden_node_for_connectivity_chromosomes_job_ACN.R $library

