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


bam_dir=`readlink -f $1`
ref=`readlink -f $2`
output_dir=`readlink -f $3`
exclude_tsv=`readlink -f $4`
mkdir -p $output_dir\/delly
mkdir -p $output_dir\/manta
declare -a pids;

tumor_bam=$bam_dir\/tumor_sorted.bam
normal_bam=$bam_dir\/normal_sorted.bam


  
LIB=`readlink -f ${BASH_SOURCE[0]} | awk '{n=split($1,f,"/"); for(i=1;i<=n-3;i++){printf "%s/", f[i]}}'`
export bam_processing_lib=$LIB\/preprocessing

cd $output_dir\/delly

$bam_processing_lib\/delly_somatic.sh $tumor_bam $normal_bam $ref 2 20 &
pids[0]=$!;

echo "delly_somatic run... delly/log/somatic"

$bam_processing_lib\/delly.sh germline $normal_bam $ref $exclude_tsv 2 20 &
pids[1]=$!;

echo "delly_germline run... delly/log/germline"

cd $output_dir\/manta
$bam_processing_lib\/manta_somatic.sh $tumor_bam $normal_bam $ref &
pids[2]=$!;

echo "manta_somatic run... manta/log"

wait ${pids[0]}
wait ${pids[1]}
wait ${pids[2]}
