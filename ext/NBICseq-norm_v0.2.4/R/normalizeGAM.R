
args <- commandArgs(TRUE)
if(length(args)!=3)
        { print("Usage: \'R --slave --args <InputData> <FragmentLength> <output> < normalize.R\'")
          q()
        }

infile = args[[1]];
fragLen = as.numeric(args[[2]]);
outfile = args[[3]];

load.success = library(mgcv,logical.return=TRUE)
if(!load.success){
	q(save="no",status=1)
	}


x = read.table(infile,header=TRUE);


fmla = "count ~ s(GC)"
nctide = names(x)[-c(1:3)]
for(nm in nctide){
	fmla = paste(fmla,"+",nm)
	}

x.model = bam(as.formula(fmla),family=poisson(),data=x)
x.model.summary = summary(x.model)
#x.model.summary$p.table

inverse_theta = max((var(x$count/x.model$fitted.values)*mean(x.model$fitted.values^2) - mean(x.model$fitted.values))/(mean(x.model$fitted.values)^2),0.1)
theta = 1/inverse_theta

x.dta.new = as.data.frame(matrix(0,nrow=2*fragLen+2,ncol=ncol(x)))
names(x.dta.new)= names(x)
x.dta.new$GC = c(0:(2*fragLen+1))/(2*fragLen+1)
row.names(x.dta.new) = paste("GC",c(0:(2*fragLen+1)),sep="")
x.dta.new.pred = predict(x.model,newdata=x.dta.new)


coef = x.model.summary$p.table
coef = as.data.frame(coef)
names(coef) = c("estimate","std","zvalue","pvalue")
coef.gc = data.frame(x.dta.new.pred,est=rep(-1,2*fragLen+2),sd=rep(-1,2*fragLen+2),pvalue=rep(x.model.summary$s.table[,4],2*fragLen+2))
names(coef.gc) = c("estimate","std","zvalue","pvalue")
coef = rbind(coef,coef.gc)

theta = c(theta,-1,-1,-1)
coef = rbind(coef,theta)
row.names(coef)[nrow(coef)] = "theta"
write.table(coef, file=outfile, row.names=TRUE,col.names=c("estimate","std","zvalue","pvalue"),quote=FALSE,sep="\t")

warnings()


