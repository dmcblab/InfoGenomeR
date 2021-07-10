#!/bin/bash
ID=$1
bin_path=$2

cat $ID | awk '{
decision=0;
if($2==1||$2==2||$2==3||$2==4||$2==5||$2==6||$2==7||$2==8||$2==9||$2==10||$2==11||$2==12||$2==13||$2==14||$2==15||$2==16||$2==17||$2==18||$2==19||$2==20||$2==21||$2==22||$2=="X"){
	cmd ="sed -n 2p '$bin_path'"$2".norm.bin"
	cmd | getline SNP_start
	close(cmd)
	split(SNP_start,start_f,sep="\t");
        cmd ="tail -n 1 '$bin_path'"$2".norm.bin"
        cmd | getline SNP_end
        close(cmd)
        split(SNP_end,end_f,sep="\t");
	if($3>start_f[1]&&$3<end_f[2])
		decision=decision+1;
}
if($4==1||$4==2||$4==3||$4==4||$4==5||$4==6||$4==7||$4==8||$4==9||$4==10||$4==11||$4==12||$4==13||$4==14||$4==15||$4==16||$4==17||$4==18||$4==19||$4==20||$4==21||$4==22||$4=="X"){
        cmd ="sed -n 2p '$bin_path'"$4".norm.bin"
        cmd | getline SNP_start
        close(cmd)
        split(SNP_start,start_f,sep="\t");
        cmd ="tail -n 1 '$bin_path'"$4".norm.bin"
        cmd | getline SNP_end
        close(cmd)
        split(SNP_end,end_f,sep="\t");
        if($5>start_f[1]&&$5<end_f[2])
		decision=decision+1;
}
if(decision==2)
	print $0 > "'$ID'.truncated"



}'
mv $ID.truncated $ID
