args <- commandArgs(TRUE)

library(ABSOLUTE)
seg.dat.fn<-read.table("copy_numbers")
norm=seg.dat.fn[!is.na(seg.dat.fn$V5) & seg.dat.fn$V5>5,];
norm=norm[norm$V6<5,]
global_copy=sum(as.numeric((2^norm$V6)*abs(norm$V4-norm$V3)))/sum(as.numeric((abs(norm$V4-norm$V3))));
global_norm=log(global_copy)/log(2);
seg.dat.fn$V6=seg.dat.fn$V6-global_norm;

seg.dat.fn<-seg.dat.fn[-1]
colnames(seg.dat.fn)<-c("Chromosome", "Start","End","Num_Probes","Segment_Mean")
levels(seg.dat.fn[,1])[levels(seg.dat.fn[,1])=="X"]<-23
seg.dat.fn[,5]=seg.dat.fn[,5];
seg.dat.fn.before=seg.dat.fn;
seg.dat.fn[,4]=abs(seg.dat.fn[,4]);
write.table(seg.dat.fn.before,"copy_numbers_ABSOLUTE_input.negative_marker", sep="\t", quote=F, col.names=T, row.names=F)

segNA=which(is.na(seg.dat.fn[,5]));
seg.dat.fn[segNA,5]=0;
write.table(seg.dat.fn,"copy_numbers_ABSOLUTE_input", sep="\t", quote=F, col.names=T, row.names=F)


sigma.p<-0.02
max.sigma.h <- 0.02
 min.ploidy <-  as.numeric(args[3]);
max.ploidy <-  as.numeric(args[4]);

if(args[5] == "T" ){
        sigma.p = 0;
        max.sigma.h = 0.01
}

########### Primary.disease important!!###
## http://software.broadinstitute.org/cancer/software/genepattern/modules/docs/ABSOLUTE##
primary.disease<-args[1]
#############
sample.name<-args[2]
results.dir<-"./ABSOLUTE_output"
max.as.seg.count <- max(15000,nrow(seg.dat.fn)+1)
max.non.clonal <- 0
 max.neg.genome <- 0
copy_num_type <- "total"

satisfied=0;
iter_ploidy_i=1;
while(satisfied < 2){
	
        optim_err=0;
	tryCatch({

		RunAbsolute("copy_numbers_ABSOLUTE_input", sigma.p, max.sigma.h, min.ploidy, max.ploidy, primary.disease, "SNP_6.0", sample.name, results.dir, max.as.seg.count, max.non.clonal, max.neg.genome, copy_num_type, maf.fn=NULL, min.mut.af=NULL, output.fn.base=NULL, verbose=TRUE)

		CreateReviewObject(args[2], paste("./ABSOLUTE_output/",args[2],".ABSOLUTE.RData",sep=""), "./ABSOLUTE_output/CreateReviewObject", "total", verbose=TRUE)



		ExtractReviewedResults(paste("./ABSOLUTE_output/CreateReviewObject/",args[2],".PP-calls_tab.txt",sep=""), "test", paste("./ABSOLUTE_output/CreateReviewObject/",args[2],".PP-modes.data.RData",sep=""), "./ABSOLUTE_output/output", args[2], "total")

		t=read.table(paste("ABSOLUTE_output/output/reviewed/SEG_MAF/",args[2],".segtab.txt",sep=""),header=T)
		for(i in 1:nrow(t)){
			if(nrow(seg.dat.fn.before[seg.dat.fn.before$Chromosome==t[i,2]&seg.dat.fn.before$Start==t[i,3]& seg.dat.fn.before$End==t[i,4],])!=0){
				t[i,5]=seg.dat.fn.before[seg.dat.fn.before$Chromosome==t[i,2]&seg.dat.fn.before$Start==t[i,3]& seg.dat.fn.before$End==t[i,4],4]
			}
			
		}
		write.table(t, paste("ABSOLUTE_output/output/reviewed/SEG_MAF/",args[2],".segtab.txt",sep=""), row.names=F, col.names=T, sep="\t", quote=F);


	}, error=function(e){optim_err<<-1;});

	satisfied=2;

	if(!exists("iter_ploidy")){
		load(paste("./ABSOLUTE_output/",args[2],".ABSOLUTE.RData",sep=""));
		temp=seg.dat$mode.res$mode.tab;
		if(nrow(temp)!=0){
			if(temp[1,"alpha"] < 0.5){
				iter_ploidy = temp[temp[,"alpha"] > 0.5,];
				for( i in 1:nrow(temp)){
					if(temp[i,"alpha"] > 0.5){
						satisfied=1;
						target_ploidy=temp[i,"tau"];
						min.ploidy = target_ploidy-0.000001;
						max.ploidy = target_ploidy+0.000001;
						iter_ploidy_i=iter_ploidy_i+1;	
						break;
					}
				}
			}
		}
	}else if(optim_err == 1){
		if(iter_ploidy_i <=nrow(iter_ploidy)){
			satisfied=1;
			target_ploidy=iter_ploidy[iter_ploidy_i,"tau"];
			min.ploidy = target_ploidy-0.000001;
			max.ploidy = target_ploidy+0.000001;
			iter_ploidy_i=iter_ploidy_i+1;
		}
	}
}

seg.dat.fn=read.table("copy_numbers_ABSOLUTE_input",header=T ,stringsAsFactors=F);
seg.dat.fn[segNA,5]=NA;
write.table(seg.dat.fn,"copy_numbers_ABSOLUTE_input", sep="\t", quote=F, col.names=T, row.names=F)

