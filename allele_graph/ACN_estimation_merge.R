new_t=data.frame();
cluster=0;
for ( i in 1:23){
        t=read.table(paste("copy_numbers.CN_opt.ACN.",i,sep=""), stringsAsFactors=F, header=T);
	lcluster=t$ACN_block;
        lcluster = unique(lcluster[complete.cases(lcluster)]);


        t$ACN_block=t$ACN_block+cluster;


	if(length(lcluster)!=0){
	        cluster=cluster+max(lcluster)+1;
	}
	new_t=rbind(new_t, t);
}

write.table(new_t,"copy_numbers.CN_opt.ACN", col.names=T, row.names=F, sep="\t", quote=F)
