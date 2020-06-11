args=commandArgs(T)
s=read.table("het_snps.format",stringsAsFactors=F)
h=read.table("hom_snps.format",stringsAsFactors=F)
names(s)[2]="coor"
names(h)[2]="coor"
#p=read.table("het_snps.format.filtered.phased_by_ACN",stringsAsFactors=F)
p=read.table("snps.format.filtered.phased_by_ACN",stringsAsFactors=F)
names(p)[3]="coor"

for(i in 1:23){
	if(i==23){
		chr="X";
	}else{
		chr=i;
	}

	m=read.table(paste(args[1],"/chr",chr,".marker",sep=""),stringsAsFactors=F)
	np=p[p[,2]==chr,];
	np=np[!duplicated(np$coor),]

	ns=s[s[,1]==chr,];
 	ns=ns[!duplicated(ns$coor),]

	nh=h[h[,1]==chr,];
 	nh=nh[!duplicated(nh$coor),]

#	nh$V4=nh$V5;

#	ns=rbind(ns,nh);
#	ns=ns[order(ns[,2]),]



#	names(ns)[2]="coor"
#	np=p[p[,2]==chr,]
# 	names(p)[3]="coor"
#	names(np)[14]="cluster"
#	names(np)[5:6] = c("phased_allele1","phased_allele2")
#
#	t=merge(ns,np, by="coor", all.x=T)
#	t= t[,c(2,1,3:9,13,14,22)];
	names(np)[14]="cluster"
	names(m)[2]="coor"
#	 t=merge(m,t, by="coor", all.x=T)
         t=merge(m,np, by="coor", all.x=T)

	t=t[,c(2,1,3,4,5,9,10,18)];

	mt = merge(t, ns, by="coor", all.x=T);
	mt[is.na(mt[,6]),6] =  mt[is.na(mt[,6]),11];
        mt[is.na(mt[,7]),7] =  mt[is.na(mt[,7]),12];

	t=mt[,c(2,1,3,4,5,6,7,8)];
	mt = merge(t, nh, by="coor", all.x=T);
        mt[is.na(mt[,6]),6] =  mt[is.na(mt[,6]),12];
        mt[is.na(mt[,7]),7] =  mt[is.na(mt[,7]),12];

        t=mt[,c(2,1,3,4,5,6,7,8)];

 

	t[is.na(t[,6]),6]=t[is.na(t[,6]),4];
        t[is.na(t[,7]),7]=t[is.na(t[,7]),4];





#	cluster=-1;
#	cluster_start=0;
#	a1=NA;
#	a2=NA;
	t$previous_allele1=NA;
	t$previous_allele2=NA;
	t$back_track=0;
	back_track=0;



#	t[!is.na(t$phased_allele1),8] = t[!is.na(t$phased_allele1),"phased_allele1"]
#        t[!is.na(t$phased_allele2),9] = t[!is.na(t$phased_allele2),"phased_allele2"]
	names(t)[6]="phased_allele1";
	names(t)[7]="phased_allele2";
	
	if(length(which(!is.na(t$cluster)))!=0){
		min_cluster = min(t$cluster, na.rm=T);
		max_cluster= max(t$cluster, na.rm=T);

		for ( cluster_i in min_cluster:max_cluster){
			w= which(!is.na(t$cluster) & t$cluster==cluster_i);
			if(length(w) > 1){
					w1=w[-1]
					w2=w[-length(w)]
					t$back_track[w1] = w1-w2;
		
					w3 = intersect ( which(t[,6] == t[,"phased_allele1"]) , w1)
					
					w4 = intersect (  which(t[,6] != t[,"phased_allele1"]) , w1)

					if(length(w3)!=0){
						t$previous_allele1[w3]=t$phased_allele1   [ w[which(w %in% w3)-1] ]
						t$previous_allele2[w3]=t$phased_allele2 [ w[which(w %in% w3)-1] ]
					}
					if(length(w4)!=0){
						t$previous_allele1[w4]= t$phased_allele2 [ w[which(w %in% w4)-1] ] 
						t$previous_allele2[w4]= t$phased_allele1 [ w[which(w %in% w4)-1] ]
					}



			}
		}
	}


#	write.table(t[,c(2,1,3:5,8:9,16:19) ],  paste("genotype.format.", chr, sep=""), row.names=F, col.names=F, sep="\t",quote=F);
        write.table(t,  paste("genotype.format.", chr, sep=""), row.names=F, col.names=F, sep="\t",quote=F);




}
