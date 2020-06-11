
#h=read.table("diploid_haplotype.11",header=T, stringsAsFactors=F)

bt=read.table("copy_numbers.CN_opt.ACN", header=T, stringsAsFactors=F);
bs=read.table("snps.format.filtered.phased_by_ACN", stringsAsFactors=F)
bs[bs[,2]=="X",2]=23;

for(chr_i in 1:23){

	if(chr_i==23){
		chr="X";
	}else{
		chr=chr_i;
	}

	h=read.table(paste("diploid_haplotype.", chr, sep=""), header=F, stringsAsFactors=F);
	names(h)=c("chr","coor","ref","hap1", "hap2")

	t=bt[bt[,1]==chr_i,];
	s=bs[bs[,2]==chr_i,];


	t$phased=0;

	sum=0;
	numbersnps=0
	i=1;

	cluster = t$ACN_block;
	cluster = unique(cluster[complete.cases(cluster)]);

	for(c_i in 1:length(cluster)){
		where = which(t$ACN_block == cluster[c_i]);
		ls=data.frame();
		hs=data.frame();
		for(w_i in 1:length(where)){
			ls = rbind(ls, s[s$V3 > t[where[w_i],2]  & s$V3 < t[where[w_i],3], ]);
			hs = rbind(hs , h[h$coor > t[where[w_i],2]  & h$coor < t[where[w_i],3] , ]);
		}
			lsh = ls[ls[,3] %in% hs$coor,];
			lsh=lsh[!duplicated(lsh[,3]),];
			hsl = hs[hs[,2] %in% ls[,3],];
			hsl=hsl[!duplicated(hsl[,2]),];
	

		
				
			if(nrow(lsh)!=0){
				 a1_hap1_cd=length(which(lsh[,5] == hsl$hap1))/nrow(lsh)
				 a1_hap2_cd=length(which(lsh[,5] == hsl$hap2))/nrow(lsh)
		
				 a2_hap1_cd=length(which(lsh[,6] == hsl$hap1))/nrow(lsh);
				 a2_hap2_cd=length(which(lsh[,6] == hsl$hap2))/nrow(lsh);
		
				if(a1_hap1_cd < a1_hap2_cd ){
					tmp=t$allele_1[where];
					t$allele_1[where]=t$allele_2[where];
					t$allele_2[where]=tmp;
					print(a1_hap2_cd);
					sum=sum+length(which(lsh[,5] == hsl$hap2));
					numbersnps=numbersnps+nrow(lsh);
				}else{
					print(a1_hap1_cd);
					sum=sum+length(which(lsh[,5] == hsl$hap1));
					numbersnps=numbersnps+nrow(lsh);


				}
				t[where,"phased"]=1;
			}


		
	}
	write.table(t, paste("copy_numbers.CN_opt.ACN.phased.", chr_i, sep="") ,sep="\t",quote=F, row.names=F, col.names=T);

}

all_t=data.frame();

for(chr_i in 1:23){
	t=read.table(paste("copy_numbers.CN_opt.ACN.phased.", chr_i, sep="") , header=T, stringsAsFactors=F);
	all_t=rbind(all_t, t);
	if(chr_i==1){
		names(all_t)=names(t);
	}
}


write.table(all_t, "copy_numbers.CN_opt.ACN.phased",  quote=F, sep="\t", row.names=F, col.names=T);


