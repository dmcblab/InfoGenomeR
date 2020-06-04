t=read.table("SVs",stringsAsFactors=F)
break_thres=100;
new_t=t[0,]
while(nrow(t)!=0){
	reverse=strsplit(t[1,6],split="to");
	where=which((t$V2 == t$V2[1] &  abs(t$V3 - t$V3[1])<break_thres & t$V4 == t$V4[1] & abs(t$V5 - t$V5[1])<break_thres & t$V6 == t$V6[1])|
		(t$V4 == t$V2[1] &  abs(t$V5 - t$V3[1])<break_thres & t$V2 == t$V4[1] & abs(t$V3 - t$V5[1])<break_thres & t$V6 == paste(reverse[[1]][2],"to",reverse[[1]][1],sep=""))
		);
	new_t[nrow(new_t)+1,]=t[1,];
	t=t[-where,]
}

write.table(new_t, "SVs", quote=F, row.names=F, col.names=F, sep="\t")

t=new_t;

coor_v=c(t[,3],t[,5])
for(i in 1:(length(coor_v))){
        coor_v=c(t[,3],t[,5])
        v=data.frame()
        v$coor=integer();
        v$row_i=integer();
        v$col_i=integer();
        v_i=1;

        v=rbind(v,c(0,0,0));
        names(v)=c("coor","row_i","col_i");
        v$coor[v_i]=coor_v[i];
        v$col_i[v_i] = if(i>nrow(t)) 5 else 3;
        v$row_i[v_i] = if(i>nrow(t)) i-nrow(t) else i;
        v_i=v_i+1;
        SIGN=1;

        explored=c();
        explored=c(explored,i);
        while(SIGN==1){
        SIGN=0;
        for(j in 1:length(coor_v)){
                if(!(j %in% explored)){
                        min_s=min(v$coor)-2;
                        max_s=min(v$coor)+nrow(v)*2;
                        if(min_s<coor_v[j] && coor_v[j]<max_s){
                ##              print(j);
                                v=rbind(v,c(0,0,0));
                                v$coor[v_i]=coor_v[j];
                                v$col_i[v_i] = if(j>nrow(t)) 5 else 3;
                                v$row_i[v_i] = if(j>nrow(t)) j-nrow(t) else j;
                                v_i=v_i+1;
                ##              v=v[order(v$coor),];
                                SIGN=1
                                explored=c(explored,j);
                        }
                }
        }

        }

        if(nrow(v)>1){
                v=v[order(v$coor),];
                print(v);
                for(k in 2:nrow(v)){
                        t[v$row_i[k],v$col_i[k]]=min(v$coor)+(k-1)*2;
                }
        }
}
write.table(t,"SVs", quote=F, sep="\t", row.names=F, col.names=F)

