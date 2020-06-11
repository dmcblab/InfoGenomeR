#!/bin/bash

##pids""
##RESULT=0;

library=$1
prefiltered=$2
snp_count_path=`readlink -f $3`
updown=7

purity=`cat ABSOLUTE_output/output/reviewed/*.test.ABSOLUTE.table.txt  | tail -n 1  | cut -f4`
purity_test=`awk 'BEGIN{if("'$purity'" > 0.97){print "over"}else{print "under"}}'`
if [[ $purity_test == "over" ]];then
	purity=0.97;
fi


for i in `seq 1 23`; do
#        Rscript $library/ACN_expand_EM_gc_amp_IQRfilter_1dis_parallel_when_CN_1_speed_up_v3.R $i $purity $snp_count_path $prefiltered $updown
        Rscript $library/ACN_estimation.R $i $purity $snp_count_path $prefiltered $updown &
	pids[${i}]=$!;
done

for pid in ${pids[*]}; do
	wait $pid;
done

Rscript $library/ACN_estimation_merge.R

