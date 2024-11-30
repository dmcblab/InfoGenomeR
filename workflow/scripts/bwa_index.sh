#!/bin/bash
LIB=`readlink -f ${BASH_SOURCE[0]} | awk '{n=split($1,f,"/"); for(i=1;i<=n-3;i++){printf "%s/", f[i]}}'`
export InfoGenomeR_lib=$LIB

input_fa=$1


mkdir -p $LIB\/humandb/custom_ref

custom_ref_dir=$LIB\/humandb/custom_ref

cat $input_fa | sed 's/>chr/>/g' > $custom_ref_dir\/ref.fa


custom_fa=$custom_ref_dir\/ref.fa

bwa-mem2 index $custom_fa
