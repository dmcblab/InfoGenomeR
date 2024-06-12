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

mkdir -p $InfoGenomeR_lib/examples/fastq

cd $InfoGenomeR_lib/examples/fastq


wget -N https://zenodo.org/records/11560492/files/normal1.fq.gz
wget -N https://zenodo.org/records/11560492/files/normal2.fq.gz
wget -N https://zenodo.org/records/11560492/files/tumor1.fq.gz
wget -N https://zenodo.org/records/11560492/files/tumor2.fq.gz
