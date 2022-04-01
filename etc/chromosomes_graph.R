args=commandArgs(T);
tag=args[1]
svf=read.table("../../SVs.CN_opt.phased")
min_chr_length=2000000

library(plotrix)
draw.half.circle=function (x, y, radius, nv = 100, border = NULL, col = NA, lty = 1,
    density = NULL, angle = 45, lwd = 1, side=upper)
{
    xylim <- par("usr")
    plotdim <- par("pin")
    ymult <- getYmult()
    angle.inc <- 2 * pi/nv
    angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
    if (length(col) < length(radius))
        col <- rep(col, length.out = length(radius))
    for (circle in 1:length(radius)) {
        xv <- cos(angles) * radius[circle] + x
        yv <- sin(angles) * radius[circle] * ymult + y

	if(side=="upper"){
		w=which(yv>y-radius/8)
	}else{
		w=which(yv<y+radius/8)
	}
	xv=xv[w]
	yv=yv[w]
#	xv <- xv[1:(length(xv)/4)]
#        yv <- yv[1:(length(yv)/4)]
        polygon(xv, yv, border = border, col = col[circle], lty = lty,
            density = density, angle = angle, lwd = lwd)
    }
    invisible(list(x = xv, y = yv))
}


t=read.table("node_keys",stringsAsFactors=F)
i=1
while(i<(nrow(t)-1)){
	if(t[i,1] ==t[i+2,1] && t[i,4] == t[i+2,4] && (t[i+1,2]+1) != t[i+2,2]){
		t[i+2,2]=t[i+1,2]+1
	}
	i=i+2
}
chr=unique(t[,1])

d_origin=t

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
	w1=min(t[t[,1]==i,2])
        w2=max(t[t[,1]==i,2])
        if(length(c(w1,w2)) !=0){
                tel=c(tel,t[  t[,1]== i & t[,2] == w1, 3]);
                tel=c(tel,t[  t[,1]== i & t[,2] == w2, 3]);
                t[  t[,1]== i & t[,2] == w1,"key"]=1

        }
}

tel=t[t[,3] %in% unique(tel),3];


node_keys=t


t=read.table(paste(Sys.getenv("InfoGenomeR_lib"),"/humandb/",Sys.getenv("Ref_version"),"_cytoBand.txt",sep=""))
t=t[t$V5=="acen",]
t=t[,c(1:3,5)]
d=data.frame()
di=1
for(i in c(1:22, "X", "Y")){
        mn=min(t[t[,1]==paste("chr",i,sep="") ,2])
        mx=max(t[t[,1]==paste("chr",i,sep="") ,3])
        d[di,1]=i
        d[di,2]=mn
        d[di,3]=mx
        di=di+1
}
cent=d;

pdf("karyotypes.pdf",width=30, height=15)
###############
max_chr_size=300000000;
ymax=max_chr_size*4;
xmax=100;
#############
#plot(c(0,100),c(-50000000,max_chr_size*4), xlab="Karyotypes", ylab="", col="white",cex.axis=1.5, cex.main=1.5, cex.lab=1.5, cex.sub=1.5, frame.plot=FALSE, axes=FALSE)



chr_color=data.frame(chrom=c(1:22,"X"),color=c("#5C575F","#774697","#3153A8","#74C1A1","#56AF5A","#FCBC00","#E40852","#A39FD1","#4979BF","#9B0E1A","#EFE83D","#E2655F","#C887BD","#ACCE5D","#8562AC","#69B70A","#206B42","#D03A01","#DCDCE8","#D70000","#80C5DE","#55429C","#F9F089"),stringsAsFactors=FALSE)



#d=read.table("node_keys",stringsAsFactors=F)

d=node_keys

i=1
while(i<(nrow(d)-1)){
        if(d[i,1] ==d[i+2,1] && d[i,4] == d[i+2,4] && (d[i+1,2]+1) != d[i+2,2]){
                d[i+2,2]=d[i+1,2]+1
        }
        i=i+2
}


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
path=scan("euler_paths.0.edges", what="character", sep="\n")


i=1;
ochr=chr
plot(c(0,100),c(0,max_chr_size*4), main="Karyotypes", xlab="",ylab="", col="white",cex.axis=1.5, cex.main=1.5, cex.lab=1.5, cex.sub=1.5, frame.plot=FALSE, axes=FALSE)



x.coor=c(0,15,30,45,60,75,0,15,30,45,60,75,0,15,30,45,60,75,0,15,30,45,60,90)
y.coor=c(rep(2.8*max_chr_size,6),rep(1.7*max_chr_size,6),rep(0.7*max_chr_size,6), rep(0,6))
chr_start=data.frame(chr=c(ochr,"der"), x.coor=x.coor[1:(length(ochr)+1)], y.coor=y.coor[1:(length(ochr)+1)])
chr_start$chr=as.character(chr_start$chr)

#chr_start=data.frame(chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","der"), x.coor=c(0,15,30,45,60,75,0,15,30,45,60,75,0,15,30,45,60,75,0,15,30,45,60,90), y.coor=c(rep(2.8*max_chr_size,6),rep(1.7*max_chr_size,6),rep(0.7*max_chr_size,6),rep(0,6)))
##x_start=1;

lchr_start=chr_start[chr_start$chr%in%ochr,]
for(lchr_i in 1:nrow(lchr_start)){
	text(lchr_start[lchr_i,2]+4, lchr_start[lchr_i,3]-50000000,paste("chromosome",lchr_start[lchr_i,1],sep=""))
}


wi=1
path_l=c()
while(wi <=length(path)){
	chr_per=data.frame(chr=c(ochr,"der"), length=rep(0,length(ochr)+1))
#        chr_per=data.frame(chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","der"), length=rep(0,24))
        f=strsplit(path[wi],split="  ")

        length_sum=0;
        j=2
        while(j<length(f[[1]])){
                nodeidx=strsplit(f[[1]][j],split="-")
                chrlength=d_origin[d_origin$node==max(nodeidx[[1]]),2]-d_origin[d_origin$node==min(nodeidx[[1]]),2]
                chr=as.character(d_origin[d_origin$node==max(nodeidx[[1]]),1]);
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


        if(length_sum<min_chr_length){
		path_l[wi]=length_sum
		wi=wi+1
                next;
        }
	wi=wi+1



	if(chr_coor=="der"){
		#print(i)
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
                       draw.half.circle(x_start+0.5,y_start,radius=0.5, col=chr_color[chr_color$chrom==d[d$node==max(nodeidx[[1]]),1],2], border=NA,nv=100000,side="lower")
                        if(as.numeric(nodeidx[[1]][1]) %% 2== 0){
				ter="pter"
			}else{
				ter="qter"
			}
#			if(j==2){
				text(x_start+0.5, y_start-25000000, ter);
#			}else{
 #                               text(x_start+0.5, y_start+25000000, ter);
#			}
                }
                if(nodeidx[[1]][2] %in% tel){
                       draw.half.circle(x_start+0.5,y_start+chrlength,radius=0.5, col=chr_color[chr_color$chrom==d[d$node==max(nodeidx[[1]]),1],2], border=NA,side="upper")
                        if(as.numeric(nodeidx[[1]][2]) %% 2== 0){
                                ter="pter"
                        }else{
                                ter="qter"
                        }
#                        if(j==2){
#                                text(x_start+0.5, y_start+chrlength-25000000,ter);
#                        }else{
                                text(x_start+0.5, y_start+chrlength+25000000, ter);
 #                       }
                }

                y_start=y_start+chrlength;

                j=j+2;
                prev_nodeidx=nodeidx;

        }
       if(!(length_sum<min_chr_length)){

                i=i+1;
                chr_start[chr_start$chr==chr_coor,2]=chr_start[chr_start$chr==chr_coor,2]+2;

	}


}


i=1;
x.coor=c(0,15,30,45,60,75,0,15,30,45,60,75,0,15,30,45,60,75,0,15,30,45,60,90)
y.coor=c(rep(2.8*max_chr_size,6),rep(1.7*max_chr_size,6),rep(0.7*max_chr_size,6), rep(0,6))
chr_start=data.frame(chr=c(ochr,"der"), x.coor=x.coor[1:(length(ochr)+1)], y.coor=y.coor[1:(length(ochr)+1)])
chr_start$chr=as.character(chr_start$chr)

#chr_start=data.frame(chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","der"), x.coor=c(0,15,30,45,60,75,0,15,30,45,60,75,0,15,30,45,60,75,0,15,30,45,60,90), y.coor=c(rep(2.8*max_chr_size,6),rep(1.7*max_chr_size,6),rep(0.7*max_chr_size,6),rep(0,6)))
##x_start=1;

wi=1
while(wi <=length(path)){
        chr_per=data.frame(chr=c(ochr,"der"), length=rep(0,length(ochr)+1))
#	chr_per=data.frame(chr=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,"X","der"), length=rep(0,24))
	f=strsplit(path[wi],split="  ")
	
	length_sum=0;
	j=2
        while(j<length(f[[1]])){
                nodeidx=strsplit(f[[1]][j],split="-")
                chrlength=d_origin[d_origin$node==max(nodeidx[[1]]),2]-d_origin[d_origin$node==min(nodeidx[[1]]),2]
                chr=as.character(d_origin[d_origin$node==max(nodeidx[[1]]),1]);
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




        rect_coor=list()
        circle_coor=list()
        circle_color=list()
        l_i=1
        seg_coor=list()
        seg_color=list()
        s_i=1
        sv_coor=list()
	sv_coor_name=list()
        v_i=1

	wi=wi+1;
	if(length_sum<min_chr_length){
#		i=i+1;
#		chr_start[chr_start$chr==chr_coor,2]=chr_start[chr_start$chr==chr_coor,2]+2;
		next;
	}

	j=2;
	x_start=chr_start[chr_start$chr==chr_coor,2];
	y_start=chr_start[chr_start$chr==chr_coor,3];
	
	while(j<length(f[[1]])){
		
		
		nodeidx=strsplit(f[[1]][j],split="-")
		
		chrlength=d[d$node==max(nodeidx[[1]]),2]-d[d$node==min(nodeidx[[1]]),2]
		chr=as.character(d[d$node==max(nodeidx[[1]]),1])
		chr_per[chr_per$chr==chr,2]=chr_per[chr_per$chr==chr,2]+chrlength;



		seg_coor[[s_i]]=c(x_start,y_start,x_start+1,y_start+chrlength)
		seg_color[[s_i]]=chr_color[chr_color$chrom==d[d$node==max(nodeidx[[1]]),1],2]
		s_i=s_i+1


                if (cent[cent[,1]==chr,3]  < d[d$node==min(nodeidx[[1]]),2] || cent[cent[,1]==chr,2] > d[d$node==max(nodeidx[[1]]),2]){


                }else{	
			#print(nodeidx[[1]]);
			cent_coor1=max(cent[cent[,1]==chr,2],   d[d$node==min(nodeidx[[1]]),2] );
			cent_coor2=min(cent[cent[,1]==chr,3],   d[d$node==max(nodeidx[[1]]),2] );
			if(nodeidx[[1]][1] < nodeidx[[1]][2]){
				cent_start= y_start+cent_coor1-d[d$node==min(nodeidx[[1]]),2];
				cent_end= cent_start + cent_coor2-cent_coor1;
			}else{
				cent_start=y_start+  d[d$node==max(nodeidx[[1]]),2] - cent_coor2;
				cent_end=cent_start+ cent_coor2-cent_coor1;
#			print(cent_start-y_start);
#			print(cent_coor2);
			}
			etag=0
			if(length(rect_coor)>0){
				for(sl_i in 1:length(rect_coor)){
					if (rect_coor[[sl_i]][2]< (cent_start+cent_end)/2 && (cent_start+cent_end)/2< rect_coor[[sl_i]][4]){
						etag=1
					}
				}

			}
			if(etag==0){
				rect_coor[[l_i]]= c(x_start, (cent_start+cent_end)/2-9000000,  x_start+1, (cent_start+cent_end)/2+9000000)
				circle_coor[[l_i]]= c(x_start+0.5,(cent_start+cent_end)/2+10000000)
				circle_color[[l_i]]=chr_color[chr_color$chrom==d[d$node==max(nodeidx[[1]]),1],2]
				l_i=l_i+1
			}
#                        rect(x_start, (cent_start+cent_end)/2-10000000, x_start+1, (cent_start+cent_end)/2+10000000, col="white", border=NA)
#                        draw.circle(x_start+0.5,(cent_start+cent_end)/2+10000000,radius=0.5, col=chr_color[chr_color$chrom==d[d$node==max(nodeidx[[1]]),1],2], border=NA,nv=100000)

#			rect(x_start, cent_start, x_start+1, cent_end, col="white", border=NA, density=20,lwd=0.8)

                }



#		grid.roundrect(x_start/xmax+0.05, (y_start+chrlength/2 )/ymax, 0.01, chrlength/ymax, gp = gpar(fill=chr_color[chr_color$chrom==d[d$node==max(nodeidx[[1]]),1],2] ), r=0);
		if(j!=2){
			is_sv=F
			is_sv_c=as.numeric(c(prev_nodeidx[[1]][1], prev_nodeidx[[1]][2], nodeidx[[1]][1] , nodeidx[[1]][2]))
			ref_c=c(is_sv_c[1]:is_sv_c[4])
#			print(is_sv_c)
#			print("h")
#			print(ref_c)
			if(length(ref_c) == 4 && all(is_sv_c == ref_c)){
				is_sv=F
			}else{
				is_sv=T
			}
			if(is_sv){
#			if(abs (as.numeric(prev_nodeidx[[1]][2]) - as.numeric(nodeidx[[1]][1])) !=1){
				sv_coor[[v_i]]=c(x_start,y_start,x_start+1, y_start)
				svf_chr1=d[d$node==nodeidx[[1]][1],1]
				coor1=d[d$node==nodeidx[[1]][1],2]
				svf_chr2=d[d$node==prev_nodeidx[[1]][2],1]
				coor2=d[d$node==prev_nodeidx[[1]][2],2]

				svf_name1=svf[(svf[,2]==svf_chr1 & svf[,3]==coor1 & svf[,4]==svf_chr2 & svf[,5]==coor2) | (svf[,4]==svf_chr1 & svf[,5]==coor1 & svf[,2]==svf_chr2 & svf[,3]==coor2),1]
#				print(svf_name1)
				sv_coor_name[[v_i]]=svf_name1
				v_i=v_i+1
			}
		}
		##print(chr_color[chr_color$chrom==d[d$node==max(nodeidx[[1]]),1],2])
		y_start=y_start+chrlength;
	
		j=j+2;
		prev_nodeidx=nodeidx;

	}
	if(length(seg_coor)>0){
	for(l in 1:length(seg_coor)){
		rect(seg_coor[[l]][1], seg_coor[[l]][2],seg_coor[[l]][3],seg_coor[[l]][4], col=seg_color[[l]],border=NA)

	}
	}
        if(length(rect_coor)>0){
        for(l in 1:length(rect_coor)){
                       rect(rect_coor[[l]][1], rect_coor[[l]][2], rect_coor[[l]][3],rect_coor[[l]][4], col="white", border=NA)	
                       draw.half.circle(circle_coor[[l]][1], circle_coor[[l]][2], radius=0.5, col=circle_color[[l]][1], border=NA,nv=100000, side="lower")
                       draw.half.circle(circle_coor[[l]][1], circle_coor[[l]][2]-20000000, radius=0.5, col=circle_color[[l]][1], border=NA,nv=100000,side="upper")
        }
        }

	if(length(sv_coor)>0){
	for(l in 1:length(sv_coor)){
		segments(sv_coor[[l]][1], sv_coor[[l]][2], sv_coor[[l]][3], sv_coor[[l]][4], lwd=2)
		if(!is.na(tag) && tag=="tag"){
			text(sv_coor[[l]][1]-1, sv_coor[[l]][2], sv_coor_name[[l]])
		}

	}
	}

       if(!(length_sum<min_chr_length)){
		i=i+1;
		chr_start[chr_start$chr==chr_coor,2]=chr_start[chr_start$chr==chr_coor,2]+2;

	}


#chr_start[chr_start$chr==chr_coor,2]=chr_start[chr_start$chr==chr_coor,2]+2;

}


#grid.roundrect(0.5,0.5,100000000/ymax,100000000/ymax, default.units="npc",gp=gpar(fill="black"))

dev.off()

