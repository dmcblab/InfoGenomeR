
args <- commandArgs(TRUE)
if(length(args)!=2)
        { print("Usage: \'R --slave --args <InputData> <output> < normalize.R\'")
          q()
        }

infile = args[[1]];
outfile = args[[2]];
#GConly = args[[3]];

library("MASS")

#if(GConly!="Y"&&GConly!="y"&&GConly!="N"&&GConly!="n"){
#	print("<GC only> must be Y/y or N/n");
#	q();
#	}

x = read.table(infile,header=TRUE);

#if(GConly=="Y"||GConly=="y"){
#	x.model = glm(count~GC,data=x,family = quasipoisson());
#	} else { x.model = glm(count~.-pos,data=x,family = quasipoisson());
#	}

x.model = glm.nb(count~.-pos,data=x)

coef = summary(x.model)$coefficients
theta = c(summary(x.model)$theta,summary(x.model)$SE.theta,-1,-1)
write.table(rbind(coef,theta), file=outfile, row.names=TRUE,col.names=c("estimate","std","zvalue","pvalue"),quote=FALSE,sep="\t")

warnings()
