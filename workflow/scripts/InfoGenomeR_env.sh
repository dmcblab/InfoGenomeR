#!/bin/bash
unset R_HOME
LIB=`readlink -f ${BASH_SOURCE[0]} | awk '{n=split($1,f,"/"); for(i=1;i<=n-3;i++){printf "%s/", f[i]}}'`
cd $LIB\/
R CMD INSTALL ext/ABSOLUTE_1.0.6.tar.gz

cd ext
tar -xvf nbicseq-norm_v0.2.4.tar.gz
cd NBICseq-norm_v0.2.4
make
cd ../

tar -xvf nbicseq-seg_v0.7.2.tar.gz
cd NBICseq-seg_v0.7.2
make
cd ../


cd $LIB
make
