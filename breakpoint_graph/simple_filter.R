sv=read.table("SVs.CN_opt",stringsAsFactors=F)
sv=sv[sv$V15!=0,];

exclude=data.frame(matrix(,ncol=3,nrow=0));
sv$V2[sv$V2=="X"]="23";
sv$V4[sv$V4=="X"]="23";

sv$is_SEG=0;
sv[sv[,1]=="<SEG>","is_SEG"]=1;
for(i in 1:nrow(sv)){
	if(sv[i,1]=="<TRA>"){
		if(as.numeric(sv[i,2])>as.numeric(sv[i,4])){
			tmp=sv[i,2];
			sv[i,2]=sv[i,4];
			sv[i,4]=tmp;
                        tmp=sv[i,3];
                        sv[i,3]=sv[i,5];
                        sv[i,5]=tmp;
			sv[i,6]=paste(strsplit(sv[i,6],"to")[[1]][2],"to",strsplit(sv[i,6],"to")[[1]][1],sep="");
		}
	}
	if(sv[i,1]=="<SEG>"){
		if(sv[i,2] != sv[i,4]){
                        if(as.numeric(sv[i,2])>as.numeric(sv[i,4])){
                                tmp=sv[i,2];
                                sv[i,2]=sv[i,4];
                                sv[i,4]=tmp;
                                tmp=sv[i,3];
                                sv[i,3]=sv[i,5];
                                sv[i,5]=tmp;
                                sv[i,6]=paste(strsplit(sv[i,6],"to")[[1]][2],"to",strsplit(sv[i,6],"to")[[1]][1],sep="");
                        }

			sv[i,1]="<TRA>";
		}else{
			if(as.numeric(sv[i,3])>as.numeric(sv[i,5])){
				tmp=sv[i,3];
				sv[i,3]=sv[i,5];
				sv[i,5]=tmp;
				sv[i,6]=paste(strsplit(sv[i,6],"to")[[1]][2],"to",strsplit(sv[i,6],"to")[[1]][1],sep="");
			}	

			if(sv[i,6]=="3to5"){
				sv[i,1]="<DEL>";
			}else if(sv[i,6] == "5to3"){
				sv[i,1]="<DUP>";
			}else{
				sv[i,1]="<INV>";
			}
		}
	}
}

sv$V2[sv$V2=="23"]="X";
sv$V4[sv$V4=="23"]="X";


cn=read.table("copy_numbers")
break_thres=1e5;

sv$class=0;

for(i in 1:nrow(sv)){
	if(sv[i,1]=="<DEL>"){
		where=which(sv[,1]=="<DUP>" & sv[,2] == sv[i,2] & (abs(sv[,3]-sv[i,3])<break_thres & sv[,3] < sv[i,3] & abs(sv[,5]-sv[i,5])<break_thres & sv[,5] < sv[i,5]) & sv[,15] == sv[i,15])
		if(length(where)==1){
			if(which(cn[,2]==sv[i,2] & cn[,3] == sv[where,3]) == which(cn[,2]==sv[i,2] & cn[,4]== sv[i,3])  &&  
		 	(which(cn[,2]==sv[i,2] & cn[,4] == sv[where,5]) +2 == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,5]) )){
				sv$class[c(i,where)]=1;
				exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[where,3];
                                exclude[nrow(exclude),3]=sv[i,3];
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[where,5];
                                exclude[nrow(exclude),3]=sv[i,5];
				next;

			}			
		}


		where=which(sv[,1]=="<DUP>" & sv[,2] == sv[i,2] & (abs(sv[,3]-sv[i,3])<break_thres & sv[,3] > sv[i,3] & abs(sv[,5]-sv[i,5])<break_thres & sv[,5] > sv[i,5]) & sv[,15] == sv[i,15])
                if(length(where)==1){
                        if(which(cn[,2]==sv[i,2] & cn[,3] == sv[where,3]) == which(cn[,2]==sv[i,2] & cn[,4]== sv[i,3]) +2  &&
                        (which(cn[,2]==sv[i,2] & cn[,4] == sv[where,5]) == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,5])  )){
                                sv$class[c(i,where)]=1;
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,3];
                                exclude[nrow(exclude),3]=sv[where,3];
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,5];
                                exclude[nrow(exclude),3]=sv[where,5];
				next;
			}
		}

                where=which(sv[,1]=="<DUP>" & sv[,2] == sv[i,2] & (abs(sv[,3]-sv[i,3])<break_thres & sv[,3] < sv[i,3] & abs(sv[,5]-sv[i,5])<break_thres & sv[,5] > sv[i,5]) & sv[,15] == sv[i,15] )
                if(length(where)==1){
                        if(which(cn[,2]==sv[i,2] & cn[,3] == sv[where,3]) == which(cn[,2]==sv[i,2] & cn[,4]== sv[i,3])  &&
                        (which(cn[,2]==sv[i,2] & cn[,4] == sv[where,5]) == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,5])  )){
                                sv$class[c(i,where)]=1;
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[where,3];
                                exclude[nrow(exclude),3]=sv[i,3];
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,5];
                                exclude[nrow(exclude),3]=sv[where,5];
                                next;
                        }
		}
	}
}


for(i in 1:nrow(sv)){
	if(sv$class[i]!=1){
	        if(sv[i,1]=="<DEL>"){
			if(abs(sv[i,3]-sv[i,5])<break_thres && which(cn[,2]==sv[i,2] & cn[,4]==sv[i,3]) +2 == which(cn[,2]==sv[i,2] & cn[,3]==sv[i,5]) ){
	                	sv$class[i]=2;
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,3];
                                exclude[nrow(exclude),3]=sv[i,5];
	                        next;
			}
		}else if(sv[i,1]=="<DUP>"){
	                if(abs(sv[i,3]-sv[i,5])<break_thres && which(cn[,2]==sv[i,2] & cn[,3]==sv[i,3]) == which(cn[,2]==sv[i,2] & cn[,4]==sv[i,5]) ){
	                        sv$class[i]=2;
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,3];
                                exclude[nrow(exclude),3]=sv[i,5];
	                        next;
	                }
		}
	}	
}


for(i in 1:nrow(sv)){
        if(sv[i,1]=="<INV>" && sv[i,6] == "5to5"){
                where=which(sv[,1]=="<INV>" & sv[,6]=="3to3" & sv[,2] == sv[i,2] & (abs(sv[,3]-sv[i,3])<break_thres & sv[,3] > sv[i,3] & abs(sv[,5]-sv[i,5])<break_thres & sv[,5] < sv[i,5]) & sv[,15] == sv[i,15])
                if(length(where)==1){
                        if(which(cn[,2]==sv[i,2] & cn[,4] == sv[where,3]) == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,3])  &&
                        (which(cn[,2]==sv[i,2] & cn[,4] == sv[where,5]) +2 == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,5] ))){
                                sv$class[c(i,where)]=1;
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,3];
                                exclude[nrow(exclude),3]=sv[where,3];
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[where,5];
                                exclude[nrow(exclude),3]=sv[i,5];
                                next;
                        }
		}
                where=which(sv[,1]=="<INV>" & sv[,6]=="3to3" & sv[,2] == sv[i,2] & (abs(sv[,3]-sv[i,3])<break_thres & sv[,3] < sv[i,3] & abs(sv[,5]-sv[i,5])<break_thres & sv[,5] > sv[i,5]) & sv[,15] == sv[i,15] )
                if(length(where)==1){
                        if(which(cn[,2]==sv[i,2] & cn[,4] == sv[where,3]) +2 == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,3])  &&
                        (which(cn[,2]==sv[i,2] & cn[,4] == sv[where,5])  == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,5] ))){
                                sv$class[c(i,where)]=1;
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[where,3];
                                exclude[nrow(exclude),3]=sv[i,3];
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,5];
                                exclude[nrow(exclude),3]=sv[where,5];
                                next;
                        }
                }

                where=which(sv[,1]=="<INV>" & sv[,6]=="3to3" & sv[,2] == sv[i,2] & (abs(sv[,3]-sv[i,3])<break_thres & sv[,3] > sv[i,3] & abs(sv[,5]-sv[i,5])<break_thres & sv[,5] > sv[i,5]) & sv[,15] == sv[i,15])
                if(length(where)==1){
                        if(which(cn[,2]==sv[i,2] & cn[,4] == sv[where,3])  == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,3])  &&
                        (which(cn[,2]==sv[i,2] & cn[,4] == sv[where,5])  == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,5] ))){
                                sv$class[c(i,where)]=1;
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,3];
                                exclude[nrow(exclude),3]=sv[where,3];
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,5];
                                exclude[nrow(exclude),3]=sv[where,5];
                                next;
                        }
                }
                where=which(sv[,1]=="<INV>" & sv[,6]=="3to3" & sv[,2] == sv[i,2] & (abs(sv[,3]-sv[i,3])<break_thres & sv[,3] < sv[i,3] & abs(sv[,5]-sv[i,5])<break_thres & sv[,5] < sv[i,5]) & sv[,15] == sv[i,15])
                if(abs(sv[i,3]-sv[i,5])<break_thres & length(where)==1){
                        if(which(cn[,2]==sv[i,2] & cn[,4] == sv[where,3]) +2 == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,3]) &&
                           which(cn[,2]==sv[i,2] & cn[,4] == sv[where,5]) == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,3]) &&
                           which(cn[,2]==sv[i,2] & cn[,4] == sv[where,5])+2 == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,5])){
                                sv$class[c(i,where)]=2;
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[where,3];
                                exclude[nrow(exclude),3]=sv[i,5];
                                next;

			}
	
		}

	}
}


for(i in 1:nrow(sv)){
        if(sv[i,1]=="<TRA>" && sv[i,6]=="3to5"){
                where=which(sv[,1]=="<TRA>" & sv[,2] == sv[i,2] & sv[,4] == sv[i,4] & sv[,6]=="5to3" &
			    abs(sv[,3]-sv[i,3])<break_thres & sv[,3] > sv[i,3] &
			    abs(sv[,5]-sv[i,5])<break_thres & sv[,5] > sv[i,5] & sv[,15] == sv[i,15]);

		if(length(where)==1){
                        if(which(cn[,2]==sv[i,2] & cn[,3] == sv[where,3])  == which(cn[,2]==sv[i,2] & cn[,4]== sv[i,3]) +2 &&
                        (which(cn[,2]==sv[i,4] & cn[,4] == sv[where,5]) == which(cn[,2]==sv[i,4] & cn[,3]== sv[i,5]) )){
                                sv$class[c(i,where)]=1;
                                exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,3];
                                exclude[nrow(exclude),3]=sv[where,3];
                                exclude[nrow(exclude)+1,1]=sv[i,4];
                                exclude[nrow(exclude),2]=sv[i,5];
                                exclude[nrow(exclude),3]=sv[where,5];
                                next;
			}
		}
	}
        if(sv[i,1]=="<TRA>" && sv[i,6]=="3to3"){
                where=which(sv[,1]=="<TRA>" & sv[,2] == sv[i,2] & sv[,4] == sv[i,4] & sv[,6]=="5to5" &
                            abs(sv[,3]-sv[i,3])<break_thres & sv[,3] > sv[i,3] &
                            abs(sv[,5]-sv[i,5])<break_thres & sv[,5] < sv[i,5] & sv[,15] == sv[i,15]);
                if(length(where)==1){
                        if(which(cn[,2]==sv[i,2] & cn[,3] == sv[where,3])  == which(cn[,2]==sv[i,2] & cn[,4]== sv[i,3]) +2 &&
                        (which(cn[,2]==sv[i,4] & cn[,3] == sv[where,5]) == which(cn[,2]==sv[i,4] & cn[,4]== sv[i,5]) )){
                                sv$class[c(i,where)]=1;
                               exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,3];
                                exclude[nrow(exclude),3]=sv[where,3];
                                exclude[nrow(exclude)+1,1]=sv[i,4];
                                exclude[nrow(exclude),2]=sv[where,5];
                                exclude[nrow(exclude),3]=sv[i,5];

                                next;
                        }
                }
	}

        if(sv[i,1]=="<TRA>" && sv[i,6]=="5to3"){
                where=which(sv[,1]=="<TRA>" & sv[,2] == sv[i,2] & sv[,4] == sv[i,4] & sv[,6]=="3to5" &
                            abs(sv[,3]-sv[i,3])<break_thres & sv[,3] > sv[i,3] &
                            abs(sv[,5]-sv[i,5])<break_thres & sv[,5] > sv[i,5] & sv[,15] == sv[i,15]);
                if(length(where)==1){
                        if(which(cn[,2]==sv[i,2] & cn[,4] == sv[where,3])  == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,3]) &&
                        (which(cn[,2]==sv[i,4] & cn[,3] == sv[where,5]) +2  == which(cn[,2]==sv[i,4] & cn[,4]== sv[i,5]) )){
                                sv$class[c(i,where)]=1;
                               exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,3];
                                exclude[nrow(exclude),3]=sv[where,3];
                                exclude[nrow(exclude)+1,1]=sv[i,4];
                                exclude[nrow(exclude),2]=sv[i,5];
                                exclude[nrow(exclude),3]=sv[where,5];

                                next;
                        }
                }

	}
        if(sv[i,1]=="<TRA>" && sv[i,6]=="5to5"){
                where=which(sv[,1]=="<TRA>" & sv[,2] == sv[i,2] & sv[,4] == sv[i,4] & sv[,6]=="3to3" &
                            abs(sv[,3]-sv[i,3])<break_thres & sv[,3] > sv[i,3] &
                            abs(sv[,5]-sv[i,5])<break_thres & sv[,5] < sv[i,5] & sv[,15] == sv[i,15]);
                if(length(where)==1){
                        if(which(cn[,2]==sv[i,2] & cn[,4] == sv[where,3])  == which(cn[,2]==sv[i,2] & cn[,3]== sv[i,3]) &&
                        (which(cn[,2]==sv[i,4] & cn[,4] == sv[where,5]) + 2 == which(cn[,2]==sv[i,4] & cn[,3]== sv[i,5]) )){
                                sv$class[c(i,where)]=1;
                               exclude[nrow(exclude)+1,1]=sv[i,2];
                                exclude[nrow(exclude),2]=sv[i,3];
                                exclude[nrow(exclude),3]=sv[where,3];
                                exclude[nrow(exclude)+1,1]=sv[i,4];
                                exclude[nrow(exclude),2]=sv[where,5];
                                exclude[nrow(exclude),3]=sv[i,5];

                                next;
                        }
                }
        }

}
sv[sv$is_SEG!=0,1]="<SEG>";



 write.table(sv[sv$class==0,1:14],"complex_SV",quote=F,sep="\t",row.names=F,col.names=F)
 write.table(sv[sv$class!=0,c(1:14,17)],"simple_SV",quote=F,sep="\t",row.names=F,col.names=F)

 exclude$X1[exclude$X1=="X"]="23";
 exclude=exclude[order(as.numeric(exclude[,1]),as.numeric(exclude[,2])),];
 exclude$X1[exclude$X1=="23"]="X";
 write.table(exclude,"exclude",quote=F,sep="\t",col.names=F,row.names=F)

