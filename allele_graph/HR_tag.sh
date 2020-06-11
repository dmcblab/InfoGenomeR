#!/bin/bash
library=$1;
ref=$2;
echo -n "" > SVs.tag
while IFS=$'\t' read -r -a line;do
		a=$((${line[2]}-100));
		b=$((${line[2]}+100));

		c=$((${line[4]}-100));
                d=$((${line[4]}+100));

		if [[ ${line[5]} == "3to3" ]] || [[ ${line[5]} == "5to5" ]];then
			perl $library/fasta.pl $ref ${line[1]} $a $b "F" > 1.fa
			perl $library/fasta.pl $ref ${line[3]} $c $d "T" > 2.fa

		else
                        perl $library/fasta.pl $ref ${line[1]} $a $b "F" > 1.fa
                        perl $library/fasta.pl $ref ${line[3]} $c $d "F" > 2.fa
		fi

		n=`blastn -query 1.fa  -subject 2.fa -dust no   -outfmt 6 | sort -nk4,4 | tail -n 1 | awk '{if($3> 90 && $4 > 100) print $0}' | wc -l`;
	
		for each in "${line[@]}";do
			printf '%s\t' $each >> SVs.tag
		done
		printf '%s\n' $n >> SVs.tag

done < SVs
