args=commandArgs(T);
thres=as.numeric(args[3]);
t=read.table("unsatisfied.list.no_telomere_ends",stringsAsFactors=F)
sv=read.table(args[1],stringsAsFactors=F)
#sv=read.table("SVs.CN_opt.filtered",stringsAsFactors=F);

t=t[t[,5]==0,];
raw_sv=read.table(args[2],stringsAsFactors=F);
one_side= if (args[4]=="T") T else F

raw_sv$ori1=NA;
raw_sv$ori2=NA;

for(i in 1:nrow(raw_sv)){
	raw_sv$ori1[i]=strsplit(raw_sv[i,6],"to")[[1]][1];
	raw_sv$ori2[i]=strsplit(raw_sv[i,6],"to")[[1]][2];
}
if(file.exists("cluster_sv")){
	if(file.info("cluster_sv")$size!=0){
		cluster_sv=read.table("cluster_sv",stringsAsFactors=F);
	}
}


sv_rbind = function(l_sv, add_sv){

	for(add_i in 1:nrow(add_sv)){
		w=which((l_sv[,2] == add_sv[add_i,2] & l_sv[,3] == add_sv[add_i,3] )|
		(l_sv[,2] == add_sv[add_i,4] & l_sv[,3] == add_sv[add_i,5])|
	       ( l_sv[,4] == add_sv[add_i,2] & l_sv[,5] == add_sv[add_i,3])|
	       ( l_sv[,4] == add_sv[add_i,4] & l_sv[,5] == add_sv[add_i,5]) );
		if(length(w)!=0){
			l_sv = l_sv[-w,];
			l_sv = rbind(l_sv, add_sv[add_i,1:14]);
		}else{
	                l_sv = rbind(l_sv, add_sv[add_i,1:14]);
		}
	}
	return(l_sv);
}






if(nrow(t)>0){
	t$ori_tag=NA;
	t[t[,3] %% 2 == 1, "ori_tag"]=5;
	t[t[,3] %% 2 == 0, "ori_tag"] = 3;

	for(i in 1:nrow(t)){
			tmp=raw_sv[(t[i,1]==raw_sv[,2] & abs(t[i,2]-raw_sv[,3]) < thres & raw_sv[,"ori1"] == t[i,"ori_tag"]) | (t[i,1] == raw_sv[,4] & abs(t[i,2] - raw_sv[,5]) < thres & raw_sv[,"ori2"] == t[i,"ori_tag"]),];

			if(nrow(tmp)!=0){
				for( j in 1:nrow(tmp)){
					if(exists("cluster_sv")){
						tmp2=cluster_sv[ (tmp[j,2] == cluster_sv[,2] & tmp[j,3] == cluster_sv[,3]) | (tmp[j,2] == cluster_sv[,4] & tmp[j,3] == cluster_sv[,5]), 15][1];
						if(!is.na(tmp2)){
							sv=sv_rbind(sv,cluster_sv[cluster_sv[,15] == tmp2,1:14]);
						}else{
							if( one_side ||   (length(which(t[,1]== tmp[j,2]  &abs( t[,2]- tmp[j,3] ) < thres & t[,"ori_tag"] == tmp[j,"ori1"] ))!=0 &
	                                                        length(which(t[,1]== tmp[j,4]  &abs( t[,2]- tmp[j,5] ) < thres & t[,"ori_tag"] == tmp[j,"ori2"] ))!=0 )
							){
								sv=sv_rbind(sv, tmp[j,1:14]);
							}
						}

					}else{
						if(one_side || length(which(t[,1]== tmp[j,2]  &abs( t[,2]- tmp[j,3] ) < thres & t[,"ori_tag"] == tmp[j,"ori1"] ))!=0 &
							length(which(t[,1]== tmp[j,4]  &abs( t[,2]- tmp[j,5] ) < thres & t[,"ori_tag"] == tmp[j,"ori2"] ))!=0
						){
							sv=sv_rbind(sv, tmp[j,1:14]);
						}
					}
				}
			}
	}
}



t=sv;
break_thres=1000;
new_t=t[0,]
while(nrow(t)!=0){
        reverse=strsplit(t[1,6],split="to");
        where=which((t$V2 == t$V2[1] &  abs(t$V3 - t$V3[1])<break_thres & t$V4 == t$V4[1] & abs(t$V5 - t$V5[1])<break_thres & t$V6 == t$V6[1])|
                (t$V4 == t$V2[1] &  abs(t$V5 - t$V3[1])<break_thres & t$V2 == t$V4[1] & abs(t$V3 - t$V5[1])<break_thres & t$V6 == paste(reverse[[1]][2],"to",reverse[[1]][1],sep=""))
                );
        new_t[nrow(new_t)+1,]=t[1,];
        t=t[-where,]
}
sv=new_t;


write.table(sv, "SVs.added", sep="\t", quote=F, col.names=F, row.names=F)
