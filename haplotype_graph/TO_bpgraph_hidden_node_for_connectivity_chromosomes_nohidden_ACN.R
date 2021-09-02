smallseg_thres=1e7;

args <- commandArgs(TRUE)#/args=c();
args_n=length(args);

chr=args[1:(args_n-1)];
#args[args_n]="filtered.format.truncated.break_adjusted.CN_opt.17";
########
ACN=read.table("copy_numbers.CN_opt.ACN.phased", header=T, stringsAsFactors=F)
ACN[ACN[,"Chromosome"] == "23",1] ="X";
t=ACN;
t=t[t$Chromosome %in% chr,]


d1=data.frame(chrom=t[,1],key=0,node=0, degree=0, stringsAsFactors=F)
dindex=1
dnode=1;
for (i in 1:nrow(t)){
	d1[dindex,1]=t[i,1]
	d1[dindex,2]=t[i,2]
	d1[dindex,3]=dnode;
	d1[dindex,4]=t[i,"allele_1"]
	dindex=dindex+1;
        dnode=dnode+1;
	d1[dindex,1]=t[i,1]
	d1[dindex,2]=t[i,3]
	d1[dindex,3]=dnode;
        d1[dindex,4]=t[i,"allele_1"]
	dindex=dindex+1;
	dnode=dnode+1;
 }

dindex=1;
d2=data.frame(chrom=t[,1],key=0,node=0, degree=0, stringsAsFactors=F)
for (i in 1:nrow(t)){
        d2[dindex,1]=t[i,1]
        d2[dindex,2]=t[i,2]
        d2[dindex,3]=dnode;
        d2[dindex,4]=t[i,"allele_2"]
        dindex=dindex+1;
        dnode=dnode+1;
        d2[dindex,1]=t[i,1]
        d2[dindex,2]=t[i,3]
        d2[dindex,3]=dnode;
        d2[dindex,4]=t[i,"allele_2"]
        dindex=dindex+1;
	dnode=dnode+1;
 }



degree1=d1;
degree2=d2;


d1$SV_CN=0;
d2$SV_CN=0;

e=data.frame(edge=as.character(),stringsAsFactors=FALSE)
e=data.frame(matrix(,nrow=0,ncol=3))
names(e)=c("edge1","edge2","type");
eindex=1


#SV=read.table("filtered.format.truncated.break_adjusted.CN_opt.AS_SV.haplotype_phased.11",stringsAsFactors=F)
if(args[args_n] != "NULL"){
SV=read.table(args[args_n],stringsAsFactors=F)
names(SV)[16:19]=c("a11","a12","a21","a22");
SV=SV[SV$V2 %in% chr | SV$V4 %in% chr,]

for(i in 1:nrow(SV)){
	if(!is.na(SV[i,"a11"])){	
		dn1=d1[d1$chrom==as.character(SV[i,2])&d1$key==SV[i,3],3]
		acn1=SV[i,"a11"];
		d1[d1$node ==dn1, "degree"]= d1[d1$node ==dn1, "degree"] - acn1;
		d1[d1$node ==dn1, "SV_CN"]= d1[d1$node ==dn1, "SV_CN"] + acn1;
	}else{
                dn1=d2[d2$chrom==as.character(SV[i,2])&d2$key==SV[i,3],3]
		acn1=SV[i,"a12"];
                d2[d2$node ==dn1, "degree"]= d2[d2$node ==dn1, "degree"] - acn1;
                d2[d2$node ==dn1, "SV_CN"]= d2[d2$node ==dn1, "SV_CN"] + acn1;

	}
	if(!is.na(SV[i,"a21"])){
		dn2=d1[d1$chrom==as.character(SV[i,4])&d1$key==SV[i,5],3]
		acn2=SV[i,"a21"];
                d1[d1$node ==dn2, "degree"]= d1[d1$node ==dn2, "degree"] - acn2;
                d1[d1$node ==dn2, "SV_CN"]= d1[d1$node ==dn2, "SV_CN"] + acn2;

	}else{
                dn2=d2[d2$chrom==as.character(SV[i,4])&d2$key==SV[i,5],3]
		acn2=SV[i,"a22"];
                d2[d2$node ==dn2, "degree"]= d2[d2$node ==dn2, "degree"] - acn2;
                d2[d2$node ==dn2, "SV_CN"]= d2[d2$node ==dn2, "SV_CN"] + acn2;
	}

	for(SV_CN in 1:min(acn1,acn2)){
		e[eindex,"edge1"]=dn1;
		e[eindex,"edge2"]=dn2;
		e[eindex,"type"]=2;
		eindex=eindex+1;
	}

}
}


r_edge_start=2;
for(i in c(1:22,"X")){
	r_edge_end=r_edge_start+nrow(d1[d1$chrom==i,])-4
	if(r_edge_start<=r_edge_end){
		edge_index=r_edge_start;
		while(edge_index<=r_edge_end){
			if(min(d1[edge_index,4],d1[edge_index+1,4])!=0){
				for(j in 1:min(d1[edge_index,4],d1[edge_index+1,4])){
					e[eindex,"edge1"]=d1[edge_index,3];
					e[eindex,"edge2"]=d1[edge_index+1,3];
					e[eindex,"type"]=1;
					eindex=eindex+1;
					d1[edge_index,4]=d1[edge_index,4]-1;
					d1[edge_index+1,4]=d1[edge_index+1,4]-1;
				}
			}
		edge_index=edge_index+2;
		}
	}
	r_edge_start=r_edge_end+4;
}

r_edge_start=2;
for(i in c(1:22,"X")){
        r_edge_end=r_edge_start+nrow(d2[d2$chrom==i,])-4
        if(r_edge_start<=r_edge_end){
                edge_index=r_edge_start;
                while(edge_index<=r_edge_end){
                        if(min(d2[edge_index,4],d2[edge_index+1,4])!=0){
                                for(j in 1:min(d2[edge_index,4],d2[edge_index+1,4])){
                                        e[eindex,"edge1"]=d2[edge_index,3];
                                        e[eindex,"edge2"]=d2[edge_index+1,3];
                                        e[eindex,"type"]=1;
                                        eindex=eindex+1;
                                        d2[edge_index,4]=d2[edge_index,4]-1;
                                        d2[edge_index+1,4]=d2[edge_index+1,4]-1;
                                }
                        }
                edge_index=edge_index+2;
                }
        }
        r_edge_start=r_edge_end+4;
}


for (i in 1:(nrow(d1)/2)){
        j=i*2-1;
        k=j+1;
	mini=min(d1[j,4],d1[k,4]);
	if(abs (d1[k,2] - d1[j,2]) < smallseg_thres){
		d1[j,4]=d1[j,4]-mini;
	        degree1[j,4]=degree1[j,4]-mini;
	        d1[k,4]=d1[k,4]-mini;
		degree1[k,4]=degree1[k,4]-mini;
	}
}

for (i in 1:(nrow(d2)/2)){
        j=i*2-1;
        k=j+1;
        mini=min(d2[j,4],d2[k,4]);
        if(abs (d2[k,2] - d2[j,2]) < smallseg_thres){
		d2[j,4]=d2[j,4]-mini;
		degree2[j,4]=degree2[j,4]-mini;
		d2[k,4]=d2[k,4]-mini;
		degree2[k,4]=degree2[k,4]-mini;
	}
}

for (i in 1:(nrow(degree1)/2)){
        j=i*2-1;
        k=j+1;
        if(degree1[j,4]!=0){
                for(CN_index in 1:degree1[j,4]){
                        e[eindex,"edge1"]=degree1[j,3]
                        e[eindex,"edge2"]=degree1[k,3]
                        e[eindex,"type"]=0;
                        eindex=eindex+1;
                }
        }
}
for (i in 1:(nrow(degree2)/2)){
        j=i*2-1;
        k=j+1;
        if(degree2[j,4]!=0){
                for(CN_index in 1:degree2[j,4]){
                        e[eindex,"edge1"]=degree2[j,3]
                        e[eindex,"edge2"]=degree2[k,3]
                        e[eindex,"type"]=0;
                        eindex=eindex+1;
                }
        }
}
############### IF the number of odd edges are not even, error message will occur######################

hidden_e=data.frame(edge=as.character(),stringsAsFactors=FALSE)
hidden_e=data.frame(matrix(,nrow=0,ncol=3))
names(hidden_e)=c("edge1","edge2","type");
hidden_node=max(d2$node)+1;

heindex=1
d_hidden= d1[d1$degree!=0,]
if(nrow(d_hidden)!=0){
	for(i in 1:(nrow(d_hidden))){
	        while(d_hidden[i,4]!=0){
	                        hidden_e[heindex,"edge1"]=d_hidden[i,3];
	                        hidden_e[heindex,"edge2"]=hidden_node;
	                        hidden_e[heindex,"type"]=3;
	                        d_hidden[i,4]=d_hidden[i,4]-1;
	                        heindex=heindex+1;
	                }
	}
}

d_hidden= d2[d2$degree!=0,]
if(nrow(d_hidden)!=0){
	for(i in 1:(nrow(d_hidden))){
	        while(d_hidden[i,4]!=0){
	                        hidden_e[heindex,"edge1"]=d_hidden[i,3];
	                        hidden_e[heindex,"edge2"]=hidden_node;
	                        hidden_e[heindex,"type"]=3;
	                        d_hidden[i,4]=d_hidden[i,4]-1;
	                        heindex=heindex+1;
	                }
	}
}
################################################


all_e=rbind(e,hidden_e);

sr_e=all_e[all_e$type<=1,]
sr_e=sr_e[order(sr_e$edge1,sr_e$edge2),]

sv_e=all_e[all_e$type==2,]
sv_e=sv_e[order(sv_e$edge1, sv_e$edge2),]

h_e=all_e[all_e$type==3,]
h_e=h_e[order(h_e$edge1, h_e$edge2),]


all_e=rbind(sr_e,sv_e,h_e);

all_e$g_index=0;
h_start=0;
for(i in 2:nrow(all_e)){
	if(all_e$type[i]==3 && h_start == 0){
		all_e$g_index[i]=all_e$g_index[i-1]+1;
		h_start=1;
	}else if(all_e$type[i]==3){
		all_e$g_index[i]=all_e$g_index[i-1];

	}else if(all_e$edge1[i]==all_e$edge1[i-1] && all_e$edge2[i]==all_e$edge2[i-1] && all_e$type[i]==all_e$type[i-1]){
		all_e$g_index[i]=all_e$g_index[i-1];
	}else{
		all_e$g_index[i]=all_e$g_index[i-1]+1;
	}

}

all_e$index=0:(nrow(all_e)-1)


all_e$degree1=0;
all_e$degree2=0;

degree=rbind(degree1,degree2);


degree[nrow(degree)+1,"node"]=hidden_node;
degree[nrow(degree),"degree"]=(length(which(all_e$edge1==hidden_node)) + length(which(all_e$edge2==hidden_node)))/2


for(i in 1:nrow(all_e)){
	all_e$degree1[i]=degree[degree$node==all_e$edge1[i],"degree"]
        all_e$degree2[i]=degree[degree$node==all_e$edge2[i],"degree"]
}

all_e$degree1=2*all_e$degree1;
all_e$degree2=2*all_e$degree2;

all_e$edge1=all_e$edge1-1;
all_e$edge2=all_e$edge2-1;


write.table(all_e[,c("edge1","edge2","degree1","degree2","type","index","g_index")],"edge_information.txt", sep="\t", quote=F, row.names=F, col.names=F)
write.table(all_e[!(all_e$edge1==(hidden_node-1) | all_e$edge2==(hidden_node-1)),c("edge1","edge2","degree1","degree2","type","index","g_index")],"edge_information.txt.nohidden", sep="\t", quote=F, row.names=F, col.names=F)


write.table(data.frame(V1=degree$node-1, V2=degree$degree, V3=degree$degree),"degrees",quote=F, sep="\t",col.names=F, row.names=F)

d1$hap=1
d2$hap=2

d=rbind(d1,d2);
d$node=d$node-1

write.table(d[,c(1:3,6)],"node_keys",quote=F ,sep="\t", col.names=T,row.names=T)

### 오일러 사이클을 위한 것. 오일러 경로를 뽑기 위해선 하나의 edge를 제거###

################################# for additional search of SVs ############
