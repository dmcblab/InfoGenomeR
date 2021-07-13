args=c();
args=commandArgs(T)
cyto=read.table(args[2],stringsAsFactors=F)
cyto[,2] = cyto[,2] - 1e6;
cyto[,3] = cyto[,3] + 1e6;
for( i in 1:nrow(cyto)){
	cyto[i,1]=strsplit(cyto[i,1],"chr")[[1]][2]
}
cyto=cyto[cyto[,5] == "acen", ]

#cyto=read.table(cyto,stringsAsFactors=F)
cn=read.table("copy_numbers.CN_opt",stringsAsFactors=F,header=T)
cn[cn[,1]==23,1]="X"
sv=read.table("SVs.CN_opt",stringsAsFactors=F)
sv=sv[sv$V15!=0,];

 sv$ori1=0
 sv$ori2=0
exclude=data.frame(matrix(,nrow=0, ncol=3));
for(i in 1:nrow(sv)){
	sv$ori1[i]=strsplit(sv[i,6],"to")[[1]][1];
        sv$ori2[i]=strsplit(sv[i,6],"to")[[1]][2];


}
if(args[1]=="T"){
	cluster_sv = data.frame(matrix(,ncol=15,nrow=0));
	names(cluster_sv)[15]="cluster"
	cluster_i=1;
}else{
	cluster_sv=  read.table("cluster_sv", stringsAsFactors=F);
	cluster_i=max(cluster_sv[,15])+1;
        names(cluster_sv)[15]="cluster"
}
##cn=cn[cn[,2]==5,]

for(i in 1:nrow(cn)){
	if(cn[i,3]-cn[i,2]>1e5){
		next;
	}
	which_SV1_f=which(sv$V2==cn[i,1]  & sv$V3==cn[i,2] & sv$ori1==5);
        which_SV1_s=which(sv$V4==cn[i,1]  & sv$V5==cn[i,2] & sv$ori2==5);
 
	if(length(which_SV1_f)==1){
		chr1=sv[which_SV1_f,4];
		coor1=sv[which_SV1_f,5];
		ori1=sv[which_SV1_f,"ori2"];
		cn1=sv[which_SV1_f,15];

	}

	if(length(which_SV1_s)==1){
                chr1=sv[which_SV1_s,2];
                coor1=sv[which_SV1_s,3];
                ori1=sv[which_SV1_s,"ori1"];
		cn1=sv[which_SV1_s,15];
	}

        which_SV2_f=which(sv$V2==cn[i,1]  & sv$V3==cn[i,3] & sv$ori1==3);
        which_SV2_s=which(sv$V4==cn[i,1]  & sv$V5==cn[i,3] & sv$ori2==3);

        if(length(which_SV2_f)==1){
                chr2=sv[which_SV2_f,4];
                coor2=sv[which_SV2_f,5];
                ori2=sv[which_SV2_f,"ori2"];
		cn2=sv[which_SV2_f,15];
        }
        if(length(which_SV2_s)==1){
                chr2=sv[which_SV2_s,2];
                coor2=sv[which_SV2_s,3];
                ori2=sv[which_SV2_s,"ori1"];
		cn2=sv[which_SV2_s,15];
        }

	
	is_cyto = length(which(cyto[,1] == cn[i,1] & cyto[,2] <= cn[i,2] & cn[i,3] <= cyto[,3])) > 0


#	if(length(c(which_SV1_f, which_SV1_s, which_SV2_f, which_SV2_s))==2 && (cn1 == cn2 ||  is_cyto )){
        if(length(c(which_SV1_f, which_SV1_s, which_SV2_f, which_SV2_s))==2){


 #              new_obj = abs(cn[i-1,"modal_cn"] - cn[i+1, "modal_cn"]);
		new_obj = abs(round(cn[i-1,"raw_expected_cn"])  - round(cn[i+1, "raw_expected_cn"]))

               obj1=cn[i,"modal_cn"]-cn1-cn[i-1,"modal_cn"];

                if(ori1 == 5){
                        where1= which(cn[,1] == chr1 & cn[,2] == coor1);
                        obj1=obj1+ cn[where1,"modal_cn"] - cn[where1-1,"modal_cn"];
			new_obj = new_obj + cn[where1,"modal_cn"] - cn[where1-1,"modal_cn"];
                }else{
                        where1= which(cn[,1] == chr1 & cn[,3] == coor1);
                       obj1=obj1+ cn[where1,"modal_cn"] - cn[where1+1,"modal_cn"];
                       new_obj= new_obj+ cn[where1,"modal_cn"]  - cn[where1+1,"modal_cn"];
                }
               obj2=cn[i,"modal_cn"]-cn2-cn[i+1,"modal_cn"];

                if(ori2 == 5){
                        where2= which(cn[,1] == chr2 & cn[,2] == coor2);
                        obj2=obj2+ cn[where2,"modal_cn"] - cn[where2-1,"modal_cn"];
                        new_obj=new_obj+ cn[where2,"modal_cn"] - cn[where2-1,"modal_cn"];

                }else{
                        where2= which(cn[,1] == chr2 & cn[,3] == coor2);
                       obj2=obj2+ cn[where2,"modal_cn"] - cn[where2+1,"modal_cn"];
                       new_obj=new_obj+ cn[where2,"modal_cn"]  - cn[where2+1,"modal_cn"];

                }
                new_obj = new_obj - 2 * min ( obj1, obj2);

		obj1 = obj1 - cn1;
		obj2 = obj2 - cn2;


		if (obj1+obj2 < new_obj){
			next;
		}


		if(args[1]=="T"){
			cluster_sv=rbind(cluster_sv,cbind(sv[c(which_SV1_f, which_SV1_s, which_SV2_f, which_SV2_s),1:14],cluster_i));
			cluster_i=cluster_i+1;
		}else{
			temp=sv[c(which_SV1_f, which_SV1_s, which_SV2_f, which_SV2_s),];
			cluster1=cluster_sv[(cluster_sv$V2==temp[1,2] & cluster_sv$V3==temp[1,3]) | (cluster_sv$V4==temp[1,2] & cluster_sv$V5==temp[1,3]) | (cluster_sv$V2==temp[1,4] & cluster_sv$V3==temp[1,5]) | (cluster_sv$V4==temp[1,4] & cluster_sv$V5==temp[1,5]),"cluster"];
			cluster2=cluster_sv[(cluster_sv$V2==temp[2,2] & cluster_sv$V3==temp[2,3]) | (cluster_sv$V4==temp[2,2] & cluster_sv$V5==temp[2,3]) | (cluster_sv$V2==temp[2,4] & cluster_sv$V3==temp[2,5]) | (cluster_sv$V4==temp[2,4] & cluster_sv$V5==temp[2,5]),"cluster"];
			if(length(cluster1)!=0 && length(cluster2) !=0){
				cluster_sv[cluster_sv$cluster==cluster2[1],"cluster"]=cluster1[1];
			}else if(length(cluster1)!=0 && length(cluster2) ==0){
				temp=temp[2,1:14];
				temp$cluster=cluster1[1];
				cluster_sv=rbind(cluster_sv, temp);
			}else if(length(cluster1)==0 && length(cluster2) !=0){
                                temp=temp[1,1:14];
                                temp$cluster=cluster2[1];
                                cluster_sv=rbind(cluster_sv, temp);
			}else{
				temp=sv[c(which_SV1_f, which_SV1_s, which_SV2_f, which_SV2_s),1:14];
				temp$cluster=cluster_i;
                                cluster_sv=rbind(cluster_sv, temp);
				cluster_i=cluster_i+1;
			}
		}
		sv=sv[-c(which_SV1_f, which_SV1_s, which_SV2_f, which_SV2_s),];
		new_i = nrow(sv)+1;
		
		
		sv[new_i,1]="<SEG>";
		sv[new_i,2]=chr1;
		sv[new_i,3]=coor1;
		sv[new_i,4]=chr2;
		sv[new_i,5]=coor2;
		sv[new_i,6]=paste(ori1,"to",ori2,sep="");
		
		exclude[nrow(exclude)+1,1]=cn[i,1];
                exclude[nrow(exclude),2]=cn[i,2];
                exclude[nrow(exclude),3]=cn[i,3];
	}
}
write.table(cluster_sv, "cluster_sv", sep="\t", quote=F, col.names=F, row.names=F)
write.table(sv[,1:14], "SVs.CN_opt.filtered.simplified", sep="\t", quote=F, col.names=F, row.names=F)
write.table(exclude, "exclude", sep="\t", quote=F, col.names=F, row.names=F)

