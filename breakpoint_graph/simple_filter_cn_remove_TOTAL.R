args=commandArgs(T);


path=args[1];


system("mkdir -p ../cn_norm_simple");

output_path="../cn_norm_simple/"
if(file.exists("exclude")){
        if(file.info("exclude")$size!=0){
                exclude=read.table("exclude",stringsAsFactors=F);
	}
}

#exclude=read.table("exclude",stringsAsFactors=F);
for(i in 1:23){
	if(i==23){
		chr="X";
	}else{
		chr=i;
	}
	t=read.table(paste(path,chr,".norm.bin",sep=""),stringsAsFactors=F,header=T);

	if(exists("exclude")){
		exclude_l=exclude[exclude[,1]==chr,];
		remove_c=c();
		if(nrow(exclude_l)!=0){
			for(j in 1:nrow(exclude_l)){
				if(length(which(t$start>= exclude_l[j,2])) && length(which(t$end<= exclude_l[j,3]))){
					remove_start=min(which(t$start>= exclude_l[j,2]));
					 remove_end=max(which(t$end<= exclude_l[j,3]));
					if(remove_start<=remove_end){
						remove_c=c(remove_c, remove_start:remove_end);	
					}
				}
			}
		}
		if(length(remove_c)>0){
			t=t[-remove_c,];
		}
	}
	write.table(t,paste(output_path,chr,".norm.bin",sep=""),quote=F, sep="\t", row.names=F, col.names=T);
}
