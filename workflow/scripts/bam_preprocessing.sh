#!/bin/bash
unset R_HOME
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

mode=$1
bam_dir=`readlink -f $2`
ref=`readlink -f $3`
output_dir=`readlink -f $4`
exclude_tsv=`readlink -f $5`
norm_map=`readlink -f $6`


if [[ $mode != "somatic" ]] && [[ $mode != "total" ]];then
        echo "mode should be somatic or total"
        exit 1
fi


mkdir -p $output_dir\/delly
mkdir -p $output_dir\/manta
mkdir -p $output_dir\/bicseq
mkdir -p $output_dir\/snp
declare -a pids;

tumor_bam=$bam_dir\/tumor_sorted.bam
normal_bam=$bam_dir\/normal_sorted.bam

  
LIB=`readlink -f ${BASH_SOURCE[0]} | awk '{n=split($1,f,"/"); for(i=1;i<=n-3;i++){printf "%s/", f[i]}}'`
export bam_processing_lib=$LIB\/preprocessing


if [[ $mode == "somatic" ]];then

	cd $output_dir\/delly

	$bam_processing_lib\/delly_somatic.sh $tumor_bam $normal_bam $ref 2 20 &>delly_somatic.log  &
	pids[0]=$!;

	echo "delly_somatic run... delly/log delly_somatic.log"

	$bam_processing_lib\/delly.sh germline $normal_bam $ref $exclude_tsv 2 20 &>delly_germline.log &
	pids[1]=$!;

	echo "delly_germline run... delly/log delly_germline.log"

	cd $output_dir\/manta
	$bam_processing_lib\/manta_somatic.sh $tumor_bam $normal_bam $ref &>manta_somatic.log &
	pids[2]=$!;

	echo "manta_somatic run... manta/log manta_somatic.log"

	cd $output_dir\/bicseq
	$bam_processing_lib\/bicseq_preprocess.sh $tumor_bam $ref $norm_map cn_norm &>bicseq_preprocess.log &
	pids[3]=$!;
	echo "bicseq_preprocess run... bicseq_preprocess.log"

	$bam_processing_lib\/bicseq_preprocess.sh $normal_bam $ref $norm_map cn_norm_germ &>bicseq_preprocess_germline.log &
	pids[4]=$!;
	echo "bicseq_preprocess_germline run... bicseq_preprocess_germline.log"

	cd $output_dir\/snp
	$bam_processing_lib\/het_SNP_detection_somatic.sh $ref $normal_bam $tumor_bam &> snp.log &
	echo "snp run... snp.log"
	pids[5]=$!;

	wait ${pids[0]}
	wait ${pids[1]}
	wait ${pids[2]}
	wait ${pids[3]}
	wait ${pids[4]}
	wait ${pids[5]}

	cd $output_dir
	cat delly/delly.format manta/manta.format > SVs
        exit 0
fi


if [[ $mode == "total" ]];then
        cd $output_dir\/delly

        $bam_processing_lib\/delly.sh tumor $tumor_bam $ref $exclude_tsv 2 20 &>delly_total.log &
        pids[0]=$!;

        echo "delly_tumor run... delly/log delly_total.log"

        cd $output_dir\/manta
        $bam_processing_lib\/manta.sh $tumor_bam $ref &>manta_total.log &
        pids[1]=$!;

        echo "manta_somatic run... manta/log manta_total.log"

        cd $output_dir\/bicseq
        $bam_processing_lib\/bicseq_preprocess.sh $tumor_bam $ref $norm_map cn_norm &>bicseq_preprocess.log &
        pids[2]=$!;
        echo "bicseq_preprocess run... bicseq_preprocess.log"

        cd $output_dir\/snp
        $bam_processing_lib\/het_SNP_detection.sh $ref $tumor_bam &> snp.log &
        echo "snp run... snp.log"
        pids[3]=$!;

        wait ${pids[0]}
        wait ${pids[1]}
        wait ${pids[2]}
        wait ${pids[3]}

        cd $output_dir
        cat delly/delly.format manta/manta.format > SVs
        exit 0

fi
