

args <- commandArgs(TRUE)

  t_col <- function(color, percent = 1, name = NULL) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (percent)*255,
               names = name)
  return(t.col)
  }




break_thres=1000;
SV=read.table("SVs.CN_opt.phased",stringsAsFactors=F)


SV$tag=1;

t<-read.table("copy_number.CN_opt.phased", header=T)





t<-cbind("Sample.1",t)
colfunc1<-colorRampPalette(c("darkolivegreen1","green4"));
colfunc2<-colorRampPalette(c("deepskyblue","dodgerblue4"));
colfunc3<-colorRampPalette(c("white","gray"));

maxt=0;
for( i in 1:23){
	if(maxt<nrow(t[t[,2]==i,])){
		maxt=nrow(t[t[,2]==i,])
	}
}	
maxt=maxt;	

pdf("ACN.pdf",width=20, height=15)
xmax=max(t[,4])
plot(c(-10000000,xmax+(maxt+10)*10000000),c(0,120), xlab="", ylab="", col="white",cex.axis=1.5, cex.main=1.5, cex.lab=1.5, cex.sub=1.5,axes=0)
j=0
chr=1
d=data.frame(chrom=0,key=0,value=0, expanded=0)
dindex=1

xlim=xmax+(maxt+10)*10000000;

colgrad=5
text(xmax+maxt*10000000 - 105000000, 122, "hap1",cex=1.5)
text(xmax+maxt*10000000 - 125000000, 122, "hap2",cex=1.5)

for( i in 1:colgrad){
	text(xmax+maxt*10000000 - 140000000, 120 - 3.5*i+1.75, i,cex=1.5)
        rect(xmax+maxt*10000000 - 110000000, 120 - 3.5*i , xmax+maxt*10000000-100000000, 120 - 3.5*(i-1), col=colfunc1(colgrad)[i], lwd=1, border=NA);
        rect(xmax+maxt*10000000 - 130000000, 120 - 3.5*i , xmax+maxt*10000000-120000000, 120 - 3.5*(i-1), col=colfunc2(colgrad)[i], lwd=1, border=NA);
}
for( i in 1:23){
	if(i==23){
		it="X"
	}else{
		it=i
	}
        text(-10000000, i*5, it, cex= 1.5);
}









curve_plot <- function(x1,x2,t1,t2,color, h){
        if( !is.na(h) && h > 0){
             a= - abs (10/(x1-x2)^2);
        }else{
             a= + abs (10/(x1-x2)^2);

        }
     b=(t1-t2)/(x1-x2) - a*(x1+x2);
     c=t1-a*x1*x1 - b*x1;
     curve(a*x*x+b*x+c, x1, x2, add=TRUE, col=color, lwd=1)

}



for (i in 1:nrow(t)){

##	if(t[i,2]!="X"){
		chrindex=as.integer(as.character(t[i,2]))*5;
##	}else{
##		chrindex=23;
##	}
	
##        if((t[i,2])!="X"){
                if(chrindex != chr ){
                        j=1;
                } else { j=j+1;
                }
#                seg=round(t[i,18]);
		segcol="gray"
	#	if(seg>14){
	#		segcol=colfunc1(14)[seg];
	#	}else{
	#		segcol="red2";
	#	}
		
		point1index=t[i,3]+(j-1)*5000000;
		point2index=t[i,4]+(j-1)*5000000;

		if(point2index-point1index<5000000){
			point1index=point1index-1500000;
			point2index=point2index+1500000;
		}

#                if(is.na(t[i,29])==T||(i>1&&t[i,1]==t[i-1,1]&&is.na(t[i-1,29])==T)||(i<nrow(t)&&t[i,1]==t[i+1,1]&&is.na(t[i+1,29])==T) || t$balanced[i] == 1){
#		if((is.na(t[i,29])==T || (t$balanced[i]  == 1) &&
#		!(i>1 && i<nrow(t) && t[i-1,1] == t[i,1] && t[i+1,1] == t[i,1] && t[i-1,"phased"]== 1 && t[i+1,"phased"]==1)) 
#		){

		if(0){
	
##		if(is.na(t[i,26])==T||t[i,28]==1||(i>1&&t[i,1]==t[i-1,1]&&is.na(t[i-1,26])==T)||(i<nrow(t)&&t[i,1]==t[i+1,1]&&is.na(t[i+1,26])==T)){
#	                points(point1index,chrindex, pch=21, bg=segcol,cex=0.7, lwd=0.5)
#			##text(point1index,chrindex-0.2,t[i,3],cex=0.5)
#	                points(point2index,chrindex, pch=21, bg=segcol,cex=0.7, lwd=0.5)
	               ## text(point2index,chrindex-0.4,t[i,4], cex=0.5)
#	                segments(point1index,chrindex,point2index,chrindex, lwd=0.5)
			rect(point1index, chrindex -0.25 , point2index, chrindex+0.25, col=segcol, lwd=1);
#                        text(point1index,chrindex,t[i,3],cex=0.25)
#                        text(point2index,chrindex,t[i,4],cex=0.25)

			if(i>1 && t[i-1,2] == t[i,2]){
				if(!pexpanded){
					segments(px, py, point1index , chrindex , col="black", lwd=1);
				}else{
					if(t[i-1,"allele_1"]!=0)
	                                        segments(px1, py1, point1index , chrindex , col="black", lwd=1);
					if(t[i-1,"allele_2"]!=0)
 	                                       segments(px2, py2, point1index , chrindex , col="black", lwd=1);
				}
			}
			pexpanded=0;
			px=point2index;
			py=chrindex;

			expanded=0;
		}else{			

                        if(i>1 && t[i-1,2] == t[i,2]){
                                if(!pexpanded){
					if(t[i,"allele_1"]!=0)
	                                        segments(px, py, point1index , chrindex+0.5 , col="black", lwd=1);
                                        if(t[i,"allele_2"]!=0)
                                        segments(px, py, point1index , chrindex -0.5, col="black", lwd=1);

                                }else{
                                        if(t[i-1,"allele_1"]!=0 && t[i,"allele_1"]!=0)
	                                        segments(px1, py1, point1index , chrindex+0.5 , col="black", lwd=1);
                                        if(t[i-1,"allele_2"]!=0 && t[i,"allele_2"]!=0)
                                        segments(px2, py2, point1index , chrindex-0.5 , col="black", lwd=1);
                                }
                        }

			
			if(t[i,5]!=0){
	                        ACN_i=0.5;
				seg=t[i,5];
				if(seg<6){
					segcol=colfunc1(6)[seg]
				}else{
					segcol="green4";
				}
				if(i>1 && t[i-1,2]==t[i,2] && t[i-1,5] == t[i,5]){
	 	                       rect(p1_point1index, chrindex +ACN_i -0.25 , point2index, chrindex+ACN_i+0.25, col=segcol, lwd=1, border=NA);
				}else{
                                       rect(point1index, chrindex +ACN_i -0.25 , point2index, chrindex+ACN_i+0.25, col=segcol, lwd=1, border=NA);
				}
				p1_point1index=point1index
#				text(point1index,chrindex+ACN_i,t[i,3],cex=0.25)
#				text(point2index,chrindex+ACN_i,t[i,4],cex=0.25)

			       px1=point2index;
				py1=chrindex+ACN_i;
	     #                   points(point1index,chrindex+ACN_i, pch=21, bg=segcol,cex=0.7, lwd=0.5)
#        	                points(point2index,chrindex+ACN_i, pch=21, bg=segcol,cex=0.7, lwd=0.5)
 #               	        segments(point1index,chrindex+ACN_i,point2index,chrindex+ACN_i, lwd=0.5)
			}
			if(t[i,6]!=0){
	                        ACN_i=-0.5;
				seg=t[i,6];
                                if(seg<6){
                                        segcol=colfunc2(6)[seg]
                                }else{
                                        segcol="dodgerblue4";
                                }
                                if(i>1 && t[i-1,2]==t[i,2] && t[i-1,6] == t[i,6]){
				       rect(p2_point1index, chrindex +ACN_i -0.25 , point2index, chrindex+ACN_i+0.25, col=segcol, lwd=1,border=NA);
				}else{
                                       rect(point1index, chrindex +ACN_i -0.25 , point2index, chrindex+ACN_i+0.25, col=segcol, lwd=1, border=NA);
				}
				p2_point1index=point1index
#                               text(point1index,chrindex+ACN_i,t[i,3],cex=0.25)
#                                text(point2index,chrindex+ACN_i,t[i,4],cex=0.25)

				px2=point2index;
				py2=chrindex+ACN_i;
#                                points(point1index,chrindex+ACN_i, pch=21, bg=segcol,cex=0.7, lwd=0.5)
 #                               points(point2index,chrindex+ACN_i, pch=21, bg=segcol,cex=0.7, lwd=0.5)
  #                              segments(point1index,chrindex+ACN_i,point2index,chrindex+ACN_i, lwd=0.5)
			}
			pexpanded=1;
			expanded=1;



		}

		d[dindex,1]=chrindex/5;
		d[dindex,2]=t[i,3];
		d[dindex,3]=point1index;
		d[dindex,4]=expanded;
		dindex=dindex+1;
		d[dindex,1]=chrindex/5;
		d[dindex,2]=t[i,4];
		d[dindex,3]=point2index;
		d[dindex,4]=expanded;
		dindex=dindex+1;


                chr=chrindex;
#        }
}


#SV=read.table(args[1], sep="\t",stringsAsFactors=F);
for(i in 1:nrow(SV)){
	if(SV[i,2]!="X"){
		SVi2=SV[i,2];
	}else{
		SVi2=23;
	}
	if(SV[i,4]!="X"){
		SVi4=SV[i,4];
	}else{
		SVi4=23;
	}

		d1=d[d$chrom==SVi2,];
		chr1=as.integer(as.character(SVi2));
		x1=d1[d1$key==SV[i,3],3];
		x1_expanded=d1[d1$key==SV[i,3],4];
#		x1_expanded=is.na(SV[i,20]);
#		print(SV[i,3]);
		d1=d[d$chrom==SVi4,];
		chr2=as.integer(as.character(SVi4));	
		x2=d1[d1$key==SV[i,5],3];
                x2_expanded=d1[d1$key==SV[i,5],4];
#		x2_expanded=is.na(SV[i,21]);

#		print(SV[i,5]);



     a=15/(x1-x2)^2


		
if(x1_expanded==0 && x2_expanded==0){
       if(SV[i,1]=="<DUP>")
		col="#E8501F"
        if(SV[i,1]=="<DEL>")
                col="#76BF72"
        if(SV[i,1]=="<INV>"){
		if(SV[i,6]=="5to5")
			col="#FFC021"
		if(SV[i,6]=="3to3")
			col="#1F9C8B"
	}
        if(SV[i,1]=="<TRA>"){
                if(SV[i,6]=="5to5")
                        col="#FFC021"
                else if (SV[i,6]=="5to3")
			col="#E8501F"
                else if (SV[i,6]=="3to5")
			col="#76BF72"
                else
			col="#1F9C8B"
	}
		

        if(SV[i,1]=="<DUP>")
                curve(-a*(x-x1)*(x-x2)+chr1*5,x1,x2,add=TRUE, col=t_col("red3", 0.1), lwd=1)
        if(SV[i,1]=="<DEL>")
                curve(-a*(x-x1)*(x-x2)+chr1*5,x1,x2,add=TRUE, col="green4", lwd=1)
        if(SV[i,1]=="<INV>")
                if(x1==x2){
                points(x1,chr1*5+0.25, pch=21, col="red",cex=3, lwd=1)
                } else{
			if(SV[i,6]=="5to5")
	                        curve(-a*(x-x1)*(x-x2)+chr1*5,x1,x2,add=TRUE, col="blue3", lwd=1)
			else
                                curve(-a*(x-x1)*(x-x2)+chr1*5,x1,x2,add=TRUE, col="yellow3", lwd=1)
                }
        if(SV[i,1]=="<TRA>"){
		if(SV[i,6]=="5to5")
	                segments(x1,chr1*5,x2,chr2*5, lwd=1, col="blue3")
		else if (SV[i,6]=="5to3")
                        segments(x1,chr1*5,x2,chr2*5, lwd=1, col="red3")
		else if (SV[i,6]=="3to5")
                        segments(x1,chr1*5,x2,chr2*5, lwd=1, col="green3")
		else
                        segments(x1,chr1*5,x2,chr2*5, lwd=1, col="yellow3")

	}

        if(SV[i,1]=="<SEG>"){
              curve_plot(x1,x2,chr1*5,chr2*5, "black", SV[i,7]);

        }


}else{


		x1_y_coor=0;
		if(SV[i,7]!=0)
#		if(!is.na(SV[i,7]))
			x1_y_coor=0.5
#		if(!is.na(SV[i,8]))
		if(SV[i,8]!=0)
			x1_y_coor=-0.5
		x2_y_coor=0;
		if(SV[i,9]!=0)
#		if(!is.na(SV[i,9]))
                       x2_y_coor=0.5
		if(SV[i,10]!=0)
#               if(!is.na(SV[i,10]))
                       x2_y_coor=-0.5

		if(x1_expanded==0)
			x1_y_coor=0;
		if(x2_expanded==0)
			x2_y_coor=0;



		if(SV[i,1]=="<SEG>"){
			if(SV[i,1]=="<TRA>" || SV[i,2] !=SV[i,4]){
                        	segments(x1,chr1*5+x1_y_coor,x2,chr2*5+x2_y_coor, lwd=1, col=t_col("purple",SV$tag[i]))
			}else{
 	                       curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col("purple",SV$tag[i]), SV[i,7]);
			}
		}else if(SV[i,6]=="5to5"){
                        if(SV[i,2] !=SV[i,4]){
                                segments(x1,chr1*5+x1_y_coor,x2,chr2*5+x2_y_coor, lwd=1, col=t_col("#FFC021",SV$tag[i]))
                        }else{
                               curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col("#FFC021",SV$tag[i]),SV[i,7]);
                        }
		}
                else if (SV[i,6]=="5to3"){
                        if(SV[i,2] !=SV[i,4]){
                                segments(x1,chr1*5+x1_y_coor,x2,chr2*5+x2_y_coor, lwd=1, col=t_col("#E8501F",SV$tag[i]))
                        }else{
                               curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col("#E8501F",SV$tag[i]),SV[i,7]);
                        }
		}
                else if (SV[i,6]=="3to5"){
                        if(SV[i,2] !=SV[i,4]){
                                segments(x1,chr1*5+x1_y_coor,x2,chr2*5+x2_y_coor, lwd=1, col=t_col("#76BF72",SV$tag[i]))
                        }else{
                               curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col("#76BF72",SV$tag[i]),SV[i,7]);
                        }
		}
                else{
                        if(SV[i,2] !=SV[i,4]){
                                segments(x1,chr1*5+x1_y_coor,x2,chr2*5+x2_y_coor, lwd=1, col=t_col("#1F9C8B",SV$tag[i]))
                        }else{
                               curve_plot(x1,x2,chr1*5+x1_y_coor,chr2*5+x2_y_coor, col=t_col("#1F9C8B",SV$tag[i]),SV[i,7]);
                        }
		}
#	}else{
#
#		if(!is.na(SV[i,16])){
#			curve_coor=0.5;
#		}else{
#			curve_coor=-0.5;
#		}
#		if(SV[i,1]=="<SEG>"){
#                         curve(-a*(x-x1)*(x-x2)+chr1*5+curve_coor,x1,x2,add=TRUE, col="black", lwd=2)
#		}else if(SV[i,1]=="<DUP>")
#               		 curve(-a*(x-x1)*(x-x2)+chr1*5+curve_coor,x1,x2,add=TRUE, col="red3", lwd=2)
#        	if(SV[i,1]=="<DEL>")
#                	curve(-a*(x-x1)*(x-x2)+chr1*5+curve_coor,x1,x2,add=TRUE, col="green4", lwd=2)
#        	if(SV[i,1]=="<INV>")
#                	if(x1==x2){
#                	points(x1,chr1*5+0.25, pch=21, col="red",cex=3, lwd=2)
#                	} else{
#                        	if(SV[i,6]=="5to5")
#                                	curve(-a*(x-x1)*(x-x2)+chr1*5+curve_coor,x1,x2,add=TRUE, col="blue3", lwd=2)
#                        	else
#                                	curve(-a*(x-x1)*(x-x2)+chr1*5+curve_coor,x1,x2,add=TRUE, col="yellow3", lwd=2)
#                	}
#





#	}


}




}



dev.off()
