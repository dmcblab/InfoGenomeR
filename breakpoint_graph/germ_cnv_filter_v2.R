args=commandArgs(T);
germ_result=args[1];
norm_bin=args[2];
tumor_bin=args[3];
chr=args[4];

	if(chr!=23){chr_index=chr}else{chr_index="X"};

	t=read.table(paste(tumor_bin,chr_index,".norm.bin",sep=""),header=T);
	
	g=read.table(germ_result);
	g=g[g$V2==chr_index,];
	g=g[complete.cases(g$V6),];
	g=g[g$V6< -0.5 | g$V6 > 0.3, ];


	remove_c=c();
        if(nrow(g)!=0){
                for(j in 1:nrow(g)){
                        if(length(which(t$start>= g[j,3])) && length(which(t$end<= g[j,4]))){
                                remove_start=min(which(t$start>= g[j,3]));
                                 remove_end=max(which(t$end<= g[j,4]));
                                if(remove_start<=remove_end){
                                        remove_c=c(remove_c, remove_start:remove_end);
                                }
                        }
                }
        }
        if(length(remove_c)>0){
                t=t[-remove_c,];
        }
        write.table(t,paste("./bicseq_norm/",chr_index,".norm.bin", sep=""),quote=F, col.names=T, row.names=F, sep="\t");




        t=read.table(paste(norm_bin,chr_index,".norm.bin",sep=""),header=T);

        g=read.table(germ_result);
        g=g[g$V2==chr_index,];
        g=g[complete.cases(g$V6),];
        g=g[g$V6< -0.5 | g$V6 > 0.3, ];


        remove_c=c();
        if(nrow(g)!=0){
                for(j in 1:nrow(g)){
                        if(length(which(t$start>= g[j,3])) && length(which(t$end<= g[j,4]))){
                                remove_start=min(which(t$start>= g[j,3]));
                                 remove_end=max(which(t$end<= g[j,4]));
                                if(remove_start<=remove_end){
                                        remove_c=c(remove_c, remove_start:remove_end);
                                }
                        }
                }
        }
        if(length(remove_c)>0){
                t=t[-remove_c,];
        }

        write.table(t,paste("./bicseq_norm_germ/",chr_index,".norm.bin",sep=""),quote=F, col.names=T, row.names=F, sep="\t");

