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


LIB=`readlink -f ${BASH_SOURCE[0]} | awk '{n=split($1,f,"/"); for(i=1;i<=n-3;i++){printf "%s/", f[i]}}'`
export InfoGenomeR_lib=$LIB

cd $InfoGenomeR_lib/humandb

idx=0;

if [[ ! -s hg19.CRG.50bp ]];then
	wget -N https://zenodo.org/records/11607705/files/hg19.CRG.50bp.tar.gz
	tar -xvf hg19.CRG.50bp.tar.gz 
fi

if [[ ! -s GRCh37.repeatmasker ]];then
	wget -N https://zenodo.org/records/11607705/files/GRCh37.repeatmasker.gz
	zcat GRCh37.repeatmasker.gz > GRCh37.repeatmasker 
fi

if [[ ! -s haplotype_1000G ]];then
        wget -N https://zenodo.org/record/5105505/files/haplotype_1000G.tar.xz && tar -xvf haplotype_1000G.tar.xz
fi


if [[ ! -s ref ]];then
        wget -N https://zenodo.org/records/11561156/files/ref.tar.gz
        tar -xvf ref.tar.gz
fi
