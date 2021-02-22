args=commandArgs(T)
thres=2000;
t=read.table(args[1],stringsAsFactors=F)

filtered=c();

for(name in c("<DEL>", "<DUP>", "<TRA>", "<INV>")){
	where=which(t$V1==name)
	if(length(where)>thres){
		PE_thres=head(sort(t$V7[where]),length(where)-2000)[length(where)-2000];
		filtered=c(filtered,which(t$V1==name & t$V7 >= PE_thres))
		
	}else{
		filtered=c(filtered,where);
	}
}

t=t[filtered,]

write.table(t,args[1],sep="\t", quote=F, col.names=F, row.names=F)
