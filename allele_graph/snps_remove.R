snp=read.table("het_snps.format",stringsAsFactors=F)

t=read.table("exclude", stringsAsFactors=F)
remove=c();
for(i in 1:nrow(t)){
	c=which(snp$V1==t[i,1] & snp$V2 < t[i,3] & snp$V2 > t[i,2]);
	remove=c(remove,c);
}

if(length(remove)>0){
	snp=snp[-remove,];
}

write.table(snp,"het_snps.format.filtered", quote=F, col.names=F, row.names=F, sep="\t")



snp = read.table("hom_snps.format",stringsAsFactors=F);

t=read.table("exclude", stringsAsFactors=F)
remove=c();
for(i in 1:nrow(t)){
        c=which(snp$V1==t[i,1] & snp$V2 < t[i,3] & snp$V2 > t[i,2]);
        remove=c(remove,c);
}

if(length(remove)>0){
        snp=snp[-remove,];
}
write.table(snp,"hom_snps.format.filtered", quote=F, col.names=F, row.names=F, sep="\t")


