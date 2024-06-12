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

wget https://zenodo.org/record/5105505/files/haplotype_1000G.tar.xz && tar -xvf haplotype_1000G.tar.xz &
pids[0]=$!;

tar -xvf hg19.CRG.50bp.tar.gz &
pids[1]=$!;



