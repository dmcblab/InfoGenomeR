#!/bin/bash
reference=hg19.fa
norm_script=./NBICseq-norm_v0.2.4/NBICseq-norm.pl
map_file=./NBICseq-norm_v0.2.4/hg19.CRG.50bp/
read_length=100
fragment_size=350
normal_bam=normal.bam


mkdir bicseq_samtools_germ
modifiedSamtools view -U BWA,bicseq_samtools_germ/,N,N $normal_bam
#mkdir bicseq_samtools_q10
#modifiedSamtools view -U BWA,bicseq_samtools_q10/,N,N -q 10 simulated_sorted.bam


mkdir cn_norm_germ
echo -e "chromName\tfaFile\tMapFile\treadPosFile\tbinFileNorm" > norm_configFile_germ;
for i in {1..23}
do
        if [ $i == 23 ]
        then
                chr="X";
        else
                chr=$i;
        fi
        echo -e "$chr\t$reference.$chr\t$map_file\hg19.50mer.CRC.chr$chr.txt\t$PWD/bicseq_samtools_germ/$chr.seq\t$PWD/cn_norm_germ/$chr.norm.bin" >> norm_configFile_germ;

done

mkdir tmp
perl $norm_script -l $read_length -s $fragment_size norm_configFile_germ ./NB_parameters_germ --tmp tmp

