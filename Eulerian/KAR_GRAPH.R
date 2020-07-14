args=commandArgs(T);
library(plotrix)
lib=args[2];
t=read.table("node_keys",stringsAsFactors=F)
chr=unique(t[,1])

tel=c();
for( i in c(1:12,16:20,"X")){
	w1=min(t[t[,1]==i,2])
	w2=max(t[t[,1]==i,2])
	if(length(c(w1,w2)) !=0){
		tel=c(tel,t[  t[,1]== i & t[,2] == w1, 3]);
		tel=c(tel,t[  t[,1]== i & t[,2] == w2, 3]);
	}
}
for( i in c(13,14,15,21,22)){
        w2=max(t[t[,1]==i,2])
        if(length(c(w2)) !=0){
                tel=c(tel,t[  t[,1]== i & t[,2] == w2, 3]);
        }
}

tel=t[t[,3] %in% unique(tel),3];


cent=read.table(paste(args[2],"/centromere",sep=""),stringsAsFactors=F)


pdf(paste(args[1],".pdf",sep=""),width=30, height=15)
###############
max_chr_size=300000000;
ymax=max_chr_size*4;
xmax=100;
#############
plot(c(0,100),c(0,max_chr_size*4), xlab="genomes", ylab="", col="white",cex.axis=1.5, cex.main=1.5, cex.lab=1.5, cex.sub=1.5)


ACN=read.table("../SNP6_level2.txt.copynumber.filtered.segmean.CN_opt.ACN.phased", header=T, stringsAsFactors=F)
ACN[ACN[,"Chromosome"] == "23",1] ="X";
t=ACN;
t=t[t$Chromosome %in% chr,]

chr_color=data.frame(chrom=unique(ACN[,1]),color=c("#5C575F","#774697","#3153A8","#74C1A1","#56AF5A","#FCBC00","#E40852","#A39FD1","#4979BF","#9B0E1A","#EFE83D","#E2655F","#C887BD","#ACCE5D","#8562AC","#69B70A","#206B42","#D03A01","#DCDCE8","#D70000","#80C5DE","#55429C","#F9F089"),stringsAsFactors=FALSE)



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




####################### BECAUSE OF NODE INDEX IN C++ PATH FINDING###############
d=rbind(d1,d2);
d$node=d$node-1
#############################################################################

e=data.frame(edge=as.character(),stringsAsFactors=FALSE)
eindex=1

for (i in 1:(nrow(d)/2)){
        j=i*2-1;
        k=j+1;
        if(d[j,4]!=0){
                for(CN_index in 1:d[j,4]){
                        e[eindex,1]=paste(d[j,3],"-",d[k,3],sep="")
                        eindex=eindex+1;
                }
        }
}
path=scan("t", what="character", sep="\n")

i=1;
chr_start=data.frame(chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","der"), x.coor=c(0,15,30,45,60,75,0,15,30,45,60,75,0,15,30,45,60,75,0,15,30,45,60,90), y.coor=c(rep(2.8*max_chr_size,6),rep(1.7*max_chr_size,6),rep(0.7*max_chr_size,6),rep(0,6)))
##x_start=1;

while(i <=length(path)){
        chr_per=data.frame(chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","der"), length=rep(0,24))
        f=strsplit(path[i],split="  ")

        length_sum=0;
        j=2
        while(j<length(f[[1]])){
                nodeidx=strsplit(f[[1]][j],split="-")
                chrlength=d[d$node==max(nodeidx[[1]]),2]-d[d$node==min(nodeidx[[1]]),2]
                chr=as.character(d[d$node==max(nodeidx[[1]]),1]);
                chr_per[chr_per$chr==chr,2]=chr_per[chr_per$chr==chr,2]+chrlength;
                j=j+2;
                length_sum=length_sum+chrlength;
              #  print(length_sum);
        }
        if(max(chr_per[,2]) > 0.5*length_sum){
##      if(max(chr_per[,2]) > 0.5*length_sum && length_sum <  250000000 ){
                chr_coor=chr_per[chr_per$length==max(chr_per[,2]),1];
        }else{
                chr_coor="der";
        }

        j=2;
        x_start=chr_start[chr_start$chr==chr_coor,2];
        y_start=chr_start[chr_start$chr==chr_coor,3];


        while(j<length(f[[1]])){


                nodeidx=strsplit(f[[1]][j],split="-")

                chrlength=d[d$node==max(nodeidx[[1]]),2]-d[d$node==min(nodeidx[[1]]),2]
                chr=as.character(d[d$node==max(nodeidx[[1]]),1])
                chr_per[chr_per$chr==chr,2]=chr_per[chr_per$chr==chr,2]+chrlength;
                if(nodeidx[[1]][1] %in% tel){
                       draw.circle(x_start+0.5,y_start,radius=0.5, col=chr_color[chr_color$chrom==d[d$node==max(nodeidx[[1]]),1],2], border=NA,nv=100000)
                        if(as.numeric(nodeidx[[1]][1]) %% 2== 0){
				ter="pter"
			}else{
				ter="qter"
			}
			if(j==2){
				text(x_start+0.5, y_start-25000000, ter);
			}else{
                                text(x_start+0.5, y_start+25000000, ter);
			}
                }
                if(nodeidx[[1]][2] %in% tel){
                       draw.circle(x_start+0.5,y_start+chrlength,radius=0.5, col=chr_color[chr_color$chrom==d[d$node==max(nodeidx[[1]]),1],2], border=NA)
                        if(as.numeric(nodeidx[[1]][2]) %% 2== 0){
                                ter="pter"
                        }else{
                                ter="qter"
                        }
                        if(j==2){
                                text(x_start+0.5, y_start+chrlength-25000000,ter);
                        }else{
                                text(x_start+0.5, y_start+chrlength+25000000, ter);
                        }
                }

                y_start=y_start+chrlength;

                j=j+2;
                prev_nodeidx=nodeidx;

        }


i=i+1;
chr_start[chr_start$chr==chr_coor,2]=chr_start[chr_start$chr==chr_coor,2]+2;

}

i=1;
chr_start=data.frame(chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","der"), x.coor=c(0,15,30,45,60,75,0,15,30,45,60,75,0,15,30,45,60,75,0,15,30,45,60,90), y.coor=c(rep(2.8*max_chr_size,6),rep(1.7*max_chr_size,6),rep(0.7*max_chr_size,6),rep(0,6)))
##x_start=1;

while(i <=length(path)){
	chr_per=data.frame(chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","der"), length=rep(0,24))
	f=strsplit(path[i],split="  ")
	
	length_sum=0;
	j=2
        while(j<length(f[[1]])){
                nodeidx=strsplit(f[[1]][j],split="-")
                chrlength=d[d$node==max(nodeidx[[1]]),2]-d[d$node==min(nodeidx[[1]]),2]
                chr=as.character(d[d$node==max(nodeidx[[1]]),1]);
                chr_per[chr_per$chr==chr,2]=chr_per[chr_per$chr==chr,2]+chrlength;
                j=j+2;
		length_sum=length_sum+chrlength;
#		print(length_sum);
        }
	if(max(chr_per[,2]) > 0.5*length_sum){
##	if(max(chr_per[,2]) > 0.5*length_sum && length_sum <  250000000 ){
		chr_coor=chr_per[chr_per$length==max(chr_per[,2]),1];
	}else{
		chr_coor="der";
	}

	j=2;
	x_start=chr_start[chr_start$chr==chr_coor,2];
	y_start=chr_start[chr_start$chr==chr_coor,3];
	
	
	while(j<length(f[[1]])){
		
		
		nodeidx=strsplit(f[[1]][j],split="-")
		
		chrlength=d[d$node==max(nodeidx[[1]]),2]-d[d$node==min(nodeidx[[1]]),2]
		chr=as.character(d[d$node==max(nodeidx[[1]]),1])
		chr_per[chr_per$chr==chr,2]=chr_per[chr_per$chr==chr,2]+chrlength;




		rect(x_start,y_start,x_start+1,y_start+chrlength, col=chr_color[chr_color$chrom==d[d$node==max(nodeidx[[1]]),1],2], border=NA)


                if (cent[cent[,1]==chr,3]  < d[d$node==min(nodeidx[[1]]),2] || cent[cent[,1]==chr,2] > d[d$node==max(nodeidx[[1]]),2]){


                }else{	
			print(nodeidx[[1]]);
			cent_coor1=max(cent[cent[,1]==chr,2],   d[d$node==min(nodeidx[[1]]),2] );
			cent_coor2=min(cent[cent[,1]==chr,3],   d[d$node==max(nodeidx[[1]]),2] );
			if(nodeidx[[1]][1] < nodeidx[[1]][2]){
				cent_start= y_start+cent_coor1-d[d$node==min(nodeidx[[1]]),2];
				cent_end= cent_start + cent_coor2-cent_coor1;
			}else{
				cent_start=y_start+  d[d$node==max(nodeidx[[1]]),2] - cent_coor2;
				cent_end=cent_start+ cent_coor2-cent_coor1;
			}
			print(cent_start-y_start);
			print(cent_coor2);
			rect(x_start, cent_start, x_start+1, cent_end, col="white", border=NA, density=20,lwd=0.8)
			

                }



#		grid.roundrect(x_start/xmax+0.05, (y_start+chrlength/2 )/ymax, 0.01, chrlength/ymax, gp = gpar(fill=chr_color[chr_color$chrom==d[d$node==max(nodeidx[[1]]),1],2] ), r=0);
		if(j!=2){
			if(abs (as.numeric(prev_nodeidx[[1]][2]) - as.numeric(nodeidx[[1]][1])) !=1){
				segments(x_start,y_start,x_start+1, y_start, lwd=1);
			}
		}
		##print(chr_color[chr_color$chrom==d[d$node==max(nodeidx[[1]]),1],2])
		y_start=y_start+chrlength;
	
		j=j+2;
		prev_nodeidx=nodeidx;

	}


i=i+1;
chr_start[chr_start$chr==chr_coor,2]=chr_start[chr_start$chr==chr_coor,2]+2;

}


#grid.roundrect(0.5,0.5,100000000/ymax,100000000/ymax, default.units="npc",gp=gpar(fill="black"))

dev.off()

