args=commandArgs(T)
t=read.table(args[2],stringsAsFactors=F)
s=read.table(args[1],stringsAsFactors=F)


v=c()
for( i in 1:nrow(s)){

	if(length(which(s[i,2]==t[,1] & s[i,3] > t[,2] & s[i,3] < t[,3])>0)){
		v=c(v,i)
	}
	if(length(which(s[i,4]==t[,1] & s[i,5] > t[,2] & s[i,5] < t[,3])>0)){
		v=c(v,i)
	}
}

v=unique(v)

if(length(v)>0){
	s=s[-v,]
}
write.table(s, args[1], quote=F, sep="\t", row.names=F, col.names=F)

