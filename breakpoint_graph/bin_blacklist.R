args=commandArgs(T)
bin=args[1]
bin_name=args[2]
blacklist=args[3]

for(chr in c(1:22,"X")){
	v=c()
	g=read.table(paste(bin,"/",chr,".norm.bin",sep=""), header=T,stringsAsFactors=F)
	t=read.table(args[3],stringsAsFactors=F)
	t=t[t[,1]==chr,]
	
	for( i in 1:nrow(t)){
		w1=which(t[i,2]<g$start & g$start<t[i,3])
		w2=which(t[i,2]<g$end & g$end<t[i,3])
		v=c(v,w1,w2)
	}
	if(length(v)>0){
		g=g[-v,]
	}
	write.table(g, paste("./",bin_name,"/",chr,".norm.bin",sep=""), sep="\t", quote=F, row.names=F, col.names=T)
}
