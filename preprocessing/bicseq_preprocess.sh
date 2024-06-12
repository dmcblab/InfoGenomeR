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
declare -a pids;

lib=`dirname $(readlink -f ${BASH_SOURCE[0]})`

reference=$2
norm_script=`which NBICseq-norm.pl`
map_file=$3
read_length=150
fragment_size=550
bam=$1

output=$4
stat=$5

samtools stats $bam > $output.bam.stats &
pids[0]=$!;

modifiedSamtools=$lib\/../ext/samtools-0.1.7a_getUnique-0.1.3/samtools


mkdir -p $output.samtools
$modifiedSamtools view -U BWA,$output.samtools/,N,N $bam &
pids[1]=$!;

wait ${pids[0]}
wait ${pids[1]}

mkdir -p $output
echo -e "chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm" > $output.configFile;
for i in {1..23}
do
        if [ $i == 23 ]
        then
                chr="X";
        else
                chr=$i;
        fi
        echo -e "$chr\t$reference.$chr\t${map_file}/hg19.50mer.CRC.chr$chr.txt\t$PWD/${output}.samtools/$chr.seq\t$PWD/${output}/$chr.norm.bin" >> $output.configFile;

done

read_length=`cat $output.bam.stats | grep ^SN | cut -f 2- | grep "average length" | awk -F "\t" '{split($2,f,"."); print f[1]}'`
fragment_size=`cat $output.bam.stats | grep ^SN | cut -f 2- | grep "insert size average" | awk -F "\t" '{split($2,f,"."); print f[1]}'`


mkdir -p $output.tmp;
perl $norm_script -l $read_length -s $fragment_size $output.configFile $output.NB_parameters --tmp $output.tmp

