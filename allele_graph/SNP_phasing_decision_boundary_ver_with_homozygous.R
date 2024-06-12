args=commandArgs(T)
prefiltered=args[2]
#SNP=read.table("SNP.format", stringsAsFactors=F)
if(prefiltered=="T"){
        SNP=read.table("het_snps.format.filtered", stringsAsFactors=F)
	SNP_hom = read.table("hom_snps.format.filtered", stringsAsFactors=F);
}else{  
        SNP=read.table("het_snps.format", stringsAsFactors=F)
	SNP_hom = read.table("hom_snps.format", stringsAsFactors=F);
}


SNP=cbind(SNP,V10=NA, V11=NA, V12=NA);
#SNP_hom = read.table("hom_snps.format.filtered", stringsAsFactors=F);
SNP_hom = cbind(SNP_hom, V10=NA, V11=NA, V12=NA);
CN_ACN=read.table("copy_numbers.CN_opt.ACN", header=T,  stringsAsFactors=F)
SNP_phased=data.frame()
#####################
s=system("cat ABSOLUTE_output/output/reviewed/*.test.ABSOLUTE.table.txt | tail -n 1", intern=T);
purity=as.numeric(strsplit(s,split="\t")[[1]][4])
max_copy=Inf;
allele_diff_for_confident=10;
zero_depth=5;
#####################
#cluster=1;

markers=list();
for( i in 1:23){
	chr=i;
	if(i==23)
		chr="X";
	markers[[i]] = read.table(paste(args[1],"/chr",chr,".marker",sep=""), stringsAsFactors=F);
}

for(CN_ACN_i in 1:nrow(CN_ACN)){

	NB_mean=CN_ACN[CN_ACN_i,"mean"];
	NB_dis1=CN_ACN[CN_ACN_i, "dis1"];
	NB_dis2=CN_ACN[CN_ACN_i, "dis2"];


	if( !is.na(CN_ACN[CN_ACN_i,"balanced"]) && ((CN_ACN[CN_ACN_i,"allele_1"] != 0 && CN_ACN[CN_ACN_i, "allele_2"]==0) || (CN_ACN[CN_ACN_i,"allele_1"] == 0 && CN_ACN[CN_ACN_i, "allele_2"]!=0) )){

                t=SNP_hom[SNP_hom$V1==CN_ACN[CN_ACN_i,1]&SNP_hom$V2>CN_ACN[CN_ACN_i,2]&SNP_hom$V2<CN_ACN[CN_ACN_i,3],];
		if(nrow(t)!=0){
			t$V11=CN_ACN[CN_ACN_i,"allele_1"];
			t$V12=CN_ACN[CN_ACN_i, "allele_2"];
			t$cluster = CN_ACN[CN_ACN_i,"ACN_block"];
			if(CN_ACN[CN_ACN_i,"allele_1"] ==0){
				t[,4] = "N";
			}else{
				tmp=t[,4];
				t[,4]=t[,5];
				t[,5]=tmp;
                               tmp=t[,6];
                                t[,6]=t[,7];
                                t[,7]=tmp;
                               tmp=t[,8];
                                t[,8]=t[,9];
                                t[,9]=tmp;
		
				t[,5]="N";
			}
		}
                background=markers[[CN_ACN[CN_ACN_i,1]]][markers[[CN_ACN[CN_ACN_i,1]]][,2] >CN_ACN[CN_ACN_i,2] &  markers[[CN_ACN[CN_ACN_i,1]]][,2] <CN_ACN[CN_ACN_i,3],];
		if(nrow(t)!=0 && nrow(background)!=0){
			background = background [!(background[,2] %in% t[,2])  ,]
		}
                if(nrow(background)!=0){
                        background$V3=".";
                        background = cbind( background, V6=NA, V7=NA, V8=NA, V9=NA, V10=NA, V11=NA, V12=NA, cluster  = CN_ACN[CN_ACN_i, "ACN_block"]);
                        background$V5 = "N";

                        if(CN_ACN[CN_ACN_i,"allele_2"] ==0){

                        }else{
                                tmp=background[,4];
                                background[,4]=background[,5];
                                background[,5]=tmp;
                        }
                }



		merge_t=rbind(t,background);
#		t[(nrow(t)+1):(nrow(t)+nrow(background)),]=background;
		if(nrow(merge_t)!=0){
			merge_t=merge_t[order(merge_t[,2]),];
			SNP_phased=rbind(SNP_phased,merge_t);
#			SNP_phased[(nrow(SNP_phased)+1):(nrow(SNP_phased)+nrow(t)),] = t;
		}

	}else if(!is.na(CN_ACN[CN_ACN_i,"mean"]) && CN_ACN[CN_ACN_i,"balanced"]!=1){
		t=SNP[SNP$V1==CN_ACN[CN_ACN_i,1]&SNP$V2>CN_ACN[CN_ACN_i,2]&SNP$V2<CN_ACN[CN_ACN_i,3],]
		
		
		if(nrow(t)!=0){	
#		for(t_i in 1:nrow(t)){
		
			if(CN_ACN[CN_ACN_i,"allele_1"]!=0 && CN_ACN[CN_ACN_i,"allele_2"]!=0){
				if(CN_ACN[CN_ACN_i,"allele_1"]<CN_ACN[CN_ACN_i,"allele_2"]){
					allele_order_1=CN_ACN[CN_ACN_i,"allele_1"];
					allele_order_2=CN_ACN[CN_ACN_i,"allele_2"];
				}else{
					allele_order_1=CN_ACN[CN_ACN_i,"allele_2"];
                                        allele_order_2=CN_ACN[CN_ACN_i,"allele_1"];
				}
				
	#			decision_bound=NB_mean*(allele_order_1*purity + 1*(1-purity));
				
	#			boundary_diff_i=Inf;
	#			while(boundary_diff_i>0 && decision_bound <= max_copy){
         #                               decision_bound=decision_bound+1;
	#				boundary_diff_i=dnbinom(decision_bound, size=NB_dis1, mu=NB_mean*(allele_order_1*purity + 1*(1-purity)))-dnbinom(decision_bound, size=NB_dis2, mu=NB_mean*(allele_order_2*purity + 1* (1-purity)));
	#			}
				
				c1=dnbinom(t[,8], size=NB_dis1, mu=NB_mean*(allele_order_1*purity + 1*(1-purity))) * dnbinom(t[,9], size=NB_dis2, mu=NB_mean*(allele_order_2*purity + 1* (1-purity)));
				c2=dnbinom(t[,9], size=NB_dis1, mu=NB_mean*(allele_order_1*purity + 1*(1-purity))) * dnbinom(t[,8], size=NB_dis2, mu=NB_mean*(allele_order_2*purity + 1* (1-purity)));
				w1= which(c1>c2  & t[,8] < t[,9]);
				w2=which(c1<=c2 & t[,8] > t[,9]);

#				w1=which(t[,8] < decision_bound  & t[,9] > decision_bound   & abs(t[,8] - t[,9]) >= allele_diff_for_confident);
#				w2=which(t[,9] < decision_bound  & t[,8] > decision_bound   & abs(t[,8] - t[,9]) >= allele_diff_for_confident);
			
				t[w1,11]=allele_order_1;
				t[w1,12]=allele_order_2;
                                t[w2,11]=allele_order_2;
                                t[w2,12]=allele_order_1;


			}

		}
	
	
		t=t[complete.cases(t[,11]),]


		if(nrow(t)!=0){

			w1=which(t[,11]==CN_ACN[CN_ACN_i,"allele_1"]);
			w2=which(t[,11]==CN_ACN[CN_ACN_i,"allele_2"]);
	
			tmp=t[w2,4];
			t[w2,4]=t[w2,5];
			t[w2,5]=tmp;
                        tmp=t[w2,6];
                        t[w2,6]=t[w2,7];
                        t[w2,7]=tmp;
                        tmp=t[w2,8];
                        t[w2,8]=t[w2,9];
                        t[w2,9]=tmp;
                        tmp=t[w2,11];
                        t[w2,11]=t[w2,12];
                        t[w2,12]=tmp;


			print(CN_ACN[CN_ACN_i,"Chromosome"]);
			print(CN_ACN[CN_ACN_i,2]);
			print(CN_ACN[CN_ACN_i,3]);
	#		print(cluster);
			print("bound");


			t=cbind(t, cluster=CN_ACN[CN_ACN_i,"ACN_block"]);
#			if(CN_ACN_i != nrow(CN_ACN))
#				if( CN_ACN[CN_ACN_i,"Chromosome"]!=CN_ACN[CN_ACN_i+1,"Chromosome"] || (is.na(CN_ACN[CN_ACN_i+1,"balanced"]) || CN_ACN[CN_ACN_i+1, "balanced"] == 1 ) )
#					cluster=cluster+1;
#	                SNP_phased[(nrow(SNP_phased)+1):(nrow(SNP_phased)+nrow(t)),] = t;
			SNP_phased=rbind(SNP_phased,t);
		}
			
	}



}

#cluster=0;
#for(i in 1:23){
#	if(i==23){
#		chr="X";
#	}else{
#		chr=i;
#	}
#	SNP_phased[SNP_phased[,2]==chr,"cluster"] = SNP_phased[SNP_phased[,2]==chr,"cluster"] + cluster;
#
#	lcluster=SNP_phased[SNP_phased[,2]==chr,"cluster"];
#	lcluster = unique(lcluster[complete.cases(lcluster)]);
#
#	cluster=cluster+max(lcluster)+1;
#
#}


#write.table(SNP_phased, "SNP_phased_by_ACN", col.names=F, quote=F, sep="\t");
write.table(SNP_phased, "snps.format.filtered.phased_by_ACN", col.names=F, quote=F, sep="\t");

