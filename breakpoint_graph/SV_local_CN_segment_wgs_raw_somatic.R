
#germ=read.table("../CNV.output_germ", header=T,stringsAsFactors=F);
#germ_break=list();
#break_thres=10000;

#for(i in 1:23){
#	if(i!=23){
#		germ_prefix=i;
#	}else{
#		germ_prefix="X"
#	}
#	germ_break[[i]]=c(germ[germ$chrom==germ_prefix,2],germ[germ$chrom==germ_prefix,3]);
#}

args <- commandArgs(TRUE)
#args[7]="/DASstorage6/leeyh/MCF7/fastq/SRR3384112_rg_sorted_dup_removed.bam"
if(args[12]=="0"){
	chr_prefix="";
}else{
	chr_prefix="chr";
}
min_seg_length=1e3;
min_log_diff=0;
min_marker=0;
bicseq_bin_size=100;
library(DNAcopy)

germ_remove<-function(x){
        j=2;
        while(j<=nrow(x)){
 #               if(x[j-1,2]==x[j,2]){
                        if(strsplit(x[j,1],"\\.")[[1]][1]!="X"){
                                chr_index=as.numeric(as.character(strsplit(x[j,1],"\\.")[[1]][1]));
                        }else{
                                chr_index=23;
                        }
                        if(length(which(abs(germ_break[[chr_index]]-x[j,2])<break_thres))){
                                x[j-1,7]=((x[j-1,3]-x[j-1,2]+1)*x[j-1,7]+(x[j,3]-x[j,2]+1)*x[j,7])/((x[j-1,3]-x[j-1,2]+1)+(x[j,3]-x[j,2]+1));
                                x[j-1,4]=x[j-1,4]+x[j,4];
                                x[j-1,3]=x[j,3];
                                x[j:(nrow(x)-1),]=x[(j+1):nrow(x),];
                                x=x[-(nrow(x)),];
                        }else{
                                j=j+1;
                        }
#                }else{
#                                j=j+1;
#                }
        }
	return(x);
}

my.merge<-function(x){
	i<-1
	j<-2
	while(j<=nrow(x)){

		if(x[i,3]-x[i,2]<min_seg_length){
			x<-x[-i,];
		}
		else if(x[j,3]-x[j,2]<min_seg_length){
			x<-x[-j,];
		}
		else{
			if(x[i,7]-x[j,7]<min_log_diff&&x[i,7]-x[j,7]>-min_log_diff){
				x[i,3]<-x[j,3];
				x[i,4]<-x[i,4]+x[j,4];
				x[i,7]<-((x[i,3]-x[i,2])*x[i,7]+(x[j,3]-x[j,2])*x[j,7])/(x[i,3]-x[i,2]+x[j,3]-x[j,2]);
				x<-x[-j,];

			}else{
				i<-i+1;
				j<-j+1;
			}

		}
	}
	return(x);
}



#table_output=data.frame()


SV=read.table(args[1], sep="\t")
bicseq_config=data.frame();
for(i in 1:23){
	if(i!=23){
		iindex=i;
	}else{
		iindex="X";
	}
	brp=data.frame(br=integer(0),ori=integer(0), stringsAsFactors=F);
	SVindex=1;
	for(j in 1:nrow(SV)){
		if(SV[j,2]==iindex){
			brp[SVindex,1]=SV[j,3];
			brp[SVindex,2]=as.integer(strsplit(as.character(SV[j,6]),"to")[[1]][1]);
			SVindex=SVindex+1;
		}
		if(SV[j,4]==iindex){
                        brp[SVindex,1]=SV[j,5];
			brp[SVindex,2]=as.integer(strsplit(as.character(SV[j,6]),"to")[[1]][2]);

			SVindex=SVindex+1;
		}
	}
	t=read.table(paste(args[3],iindex,".norm.bin",sep=""), header=T);
	normal=read.table(paste(args[8],iindex,".norm.bin",sep=""), header=T);
	brp=rbind(c(t[1,1],5),brp,c(t[nrow(t),2],3));
        brp=brp[order(brp[1],-brp[2]),]
	#FINAL=data.frame(ID=character(), chrom=character(), loc.start=numeric(), loc.end=numeric(), num.mark=numeric(), seg.mean=numeric(), bstat=numeric(), pval=numeric(), lcl=numeric(), ucl=numeric(), "brp[l - 1, 1]"=numeric(), "brp[l, 1]"=numeric(), stringsAsFactors=F)
#	bicseq_config=data.frame();
        for(l in 2:nrow(brp)){
		tindex1=0;nindex1=0;
		tindex2=0;nindex2=0;
		tindex_t=which(brp[l-1,1]<=t[,1] & t[,2]<=brp[l,1]);
                nindex_t=which(brp[l-1,1]<=normal[,1] & normal[,2]<=brp[l,1]);
		if(length(tindex_t)!=0){
			tindex1=min(tindex_t);
			tindex2=max(tindex_t);
		}
                if(length(nindex_t)!=0){
                        nindex1=min(nindex_t);
                        nindex2=max(nindex_t);
                }

                if(brp[l-1,2]==5){
			brp_adjusted1=brp[l-1,1]
		}else{
                        brp_adjusted1=brp[l-1,1]+1
		}
                if(brp[l,2]==5){
                        brp_adjusted2=brp[l,1]-1
		}else{
                        brp_adjusted2=brp[l,1]
		}

		if(brp_adjusted1<=brp_adjusted2){
			if(tindex1<tindex2 && nindex1<nindex2){
				cn_table=t[tindex1:tindex2,]
				write.table(cn_table,paste(args[9],iindex,".norm.bin.",brp_adjusted1,"-",brp_adjusted2,sep=""),quote=F, sep="\t", col.names=T, row.names=F)
                                cn_table=normal[nindex1:nindex2,]
                                write.table(cn_table,paste(args[10],iindex,".norm.bin.",brp_adjusted1,"-",brp_adjusted2,sep=""),quote=F, sep="\t", col.names=T, row.names=F)
				bicseq_config=rbind(bicseq_config,data.frame(chromName=paste(iindex,".",brp_adjusted1,"-",brp_adjusted2,sep=""),binFileNorm.Case=paste(args[9],iindex,".norm.bin.",brp_adjusted1,"-",brp_adjusted2,sep=""),binFileNorm.Control=paste(args[10],iindex,".norm.bin.",brp_adjusted1,"-",brp_adjusted2,sep="")))
			}
		}
	}
}
	
write.table(bicseq_config, "configFile", quote=F, sep="\t", col.names=T, row.names=F);
#system(paste("perl",args[5],"--lambda",args[6],"--tmp",args[4],"--strict","--fig CNV.png","configFile","CNV.output","--noscale",sep=" "));
#system(paste("perl",args[5],"--tmp",args[4],"--fig CNV.png","configFile","CNV.output","--nrm","--noscale","--control",sep=" "));
system(paste("perl",args[5],"--tmp",args[4], "--lambda",args[6],"--fig CNV.png","configFile","CNV.output","--nrm","--noscale","--control",sep=" "));

system("cat CNV.output | awk '{print $1\"\t\"$2\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6\"\t\"$9}' > CNV.output.somatic_format");
#system(paste("perl",args[5],"--lambda",args[6],"--tmp",args[4],"--strict","--fig CNV.png","configFile","CNV.output",sep=" "));
segs2=read.table("CNV.output.somatic_format", header=T, stringsAsFactors=F);
names(segs2)=c("X1","X2","X3","X4","X5","X6","X7");
FINAL=segs2[-(1:nrow(segs2)),]
for(i in 1:23){
        if(i!=23){
                iindex=i;
        }else{
                iindex="X";
        }
        brp=data.frame(br=integer(0),ori=integer(0), stringsAsFactors=F);
        SVindex=1;
        for(j in 1:nrow(SV)){
                if(SV[j,2]==iindex){
                        brp[SVindex,1]=SV[j,3];
                        brp[SVindex,2]=as.integer(strsplit(as.character(SV[j,6]),"to")[[1]][1]);
                        SVindex=SVindex+1;
                }
                if(SV[j,4]==iindex){
                        brp[SVindex,1]=SV[j,5];
                        brp[SVindex,2]=as.integer(strsplit(as.character(SV[j,6]),"to")[[1]][2]);

                        SVindex=SVindex+1;
                }
        }
        t=read.table(paste(args[3],iindex,".norm.bin",sep=""), header=T);
        normal=read.table(paste(args[8],iindex,".norm.bin",sep=""), header=T);

        brp=rbind(c(t[1,1],5),brp,c(t[nrow(t),2],3));
        brp=brp[order(brp[1],-brp[2]),]

        for(l in 2:nrow(brp)){
                tindex1=0;nindex1=0;
                tindex2=0;nindex2=0;
                tindex_t=which(brp[l-1,1]<=t[,1] & t[,2]<=brp[l,1]);
                nindex_t=which(brp[l-1,1]<=normal[,1] & normal[,2]<=brp[l,1]);
                if(length(tindex_t)!=0){
                        tindex1=min(tindex_t);
                        tindex2=max(tindex_t);
                }
                if(length(nindex_t)!=0){
                        nindex1=min(nindex_t);
                        nindex2=max(nindex_t);
                }

                if(brp[l-1,2]==5){
                        brp_adjusted1=brp[l-1,1]
                }else{
                        brp_adjusted1=brp[l-1,1]+1
                }
                if(brp[l,2]==5){
                        brp_adjusted2=brp[l,1]-1
                }else{
                        brp_adjusted2=brp[l,1]
                }

                if(brp_adjusted1<=brp_adjusted2){
                if(tindex1<tindex2 && nindex1<nindex2){
			segs2_i=which(segs2[,1]==paste(iindex,".",brp_adjusted1,"-",brp_adjusted2,sep=""));
			if(length(segs2_i)!=0){
                                segs2_u=segs2[min(segs2_i):max(segs2_i),];
                                segs2_u=my.merge(segs2_u);
		#		segs2_u=germ_remove(segs2_u);
				segs2_u[1,2]<-brp_adjusted1;
	                        segs2_u[nrow(segs2_u),3]<-brp_adjusted2;
				print(segs2_u);
				FINAL=rbind(FINAL,segs2_u);
		
			}else{
                        	NAframe=data.frame(matrix(c(as.character(iindex),brp_adjusted1,brp_adjusted2,rep(NA,times=4)),ncol=7), stringsAsFactors=F);
                        	FINAL=rbind(FINAL,NAframe)
			}
                } else {
                        NAframe=data.frame(matrix(c(as.character(iindex),brp_adjusted1,brp_adjusted2,rep(NA,times=4)),ncol=7), stringsAsFactors=F);
			FINAL=rbind(FINAL,NAframe)
		}
		}
	}
}
for(FINAL_i in 1:nrow(FINAL)){
	FINAL[FINAL_i,1]=strsplit(FINAL[FINAL_i,1],"\\.")[[1]][1];
}
#print("debug1");
FINAL=data.frame(ID="Sample.1", chrom=FINAL[,1], loc.start=FINAL[,2], loc.end=FINAL[,3], num.mark=FINAL[,4], seg.mean=FINAL[,7], bstat=NA, pval=NA, lcl=NA, ucl=NA,"brp[l - 1, 1]"=FINAL[,2],"brp[l, 1]"=FINAL[,3])
FINAL=FINAL[order(FINAL[,2],FINAL[,3]),];
#print("debug2");
names(FINAL)<-c("ID","chrom","loc.start","loc.end","num.mark","seg.mean","bstat","pval","lcl","ucl","brp[l - 1, 1]","brp[l, 1]")

############# marker <15#################
for(FINAL_i in 1:nrow(FINAL)){
	if(is.na(FINAL[FINAL_i,5])!=T&&as.integer(FINAL[FINAL_i,5])<min_marker){
        	FINAL[FINAL_i,6]="NA"
		FINAL[FINAL_i,5]="NA"
	}
}
write.table(FINAL,"copy_numbers",sep="\t", quote=F, row.names=F, col.names=F)
#table_output<-rbind(table_output,FINAL);

#write.table(table_output, "SNP6_level2.txt.copynumber.filtered.segmean", sep="\t", quote=F, row.names=F, col.names=F)
t=read.table("copy_numbers")
t$V2=factor(t$V2, levels=c(1:22,"X"))
t=t[order(t$V2,t$V3),]
write.table(t,"copy_numbers",sep="\t", quote=F, row.names=F, col.names=F)



segmean_output=read.table("copy_numbers");
cnv_output=read.table("CNV.output.somatic_format", header=T,stringsAsFactors=F)
#global_norm=cnv_output[1,7]-log(cnv_output[1,5]/cnv_output[1,6],base=2);
new_output=segmean_output[0,];

for(j in 1:23){
        if(j!=23){
                iindex=j;
        }else{
                iindex="X";
        }	
        t=read.table(paste(args[3],iindex,".norm.bin",sep=""), header=T);
        normal=read.table(paste(args[8],iindex,".norm.bin",sep=""), header=T);

	S=segmean_output[segmean_output$V2==iindex,]
#	S$V13=0;
	for(i in 1:nrow(S)){
		if(is.na(S[i,6])){
	#		a=t[(S$V3[i]<t$start & t$start<S$V4[i]) | (S$V3[i]< t$end & t$end<S$V4[i]) | (t$start < S$V3[i] & S$V4[i] < t$end ),]
	#		if(nrow(a)!=0){
				S[i,5]=-ceiling((S$V4[i]-S$V3[i])/bicseq_bin_size);
				S[i,6]=NA;
	
		}
	}
	new_output=rbind(new_output, S);
}


write.table(new_output,"copy_numbers",sep="\t", quote=F, row.names=F, col.names=F)


