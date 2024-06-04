


infile = "tmp/binfiltering_MKbt1d8.txt"

infile = "tmp/binfiltering_soeoGeo.txt"

infile = "tmp/binfiltering__QD2XUM.txt"

if(use.loess!=TRUE){
	load.success = library(mgcv,logical.return=TRUE)
	if(!load.success){
		q(save="no",status=1)
		}
	}

load.success = library(mclust,logical.return=TRUE)
if(!load.success){
		q(save="no",status=1)
		}


cutoff = 0.5
num.sd = 3 ## if use gam, use num.sd sds to  mark the outliers
epsilon = 0.01
x_in = read.table(infile,header=T,sep="\t")

ind = c(1:nrow(x_in))
n.sample = 10000
if(nrow(x_in)>n.sample){
	gc.qt = quantile(x_in$gc,probs=c(0.001,0.999))
	ind1 = which(x_in$gc<gc.qt[1] | x_in$gc > gc.qt[2])
	ind2 = sample(ind,n.sample)
	ind3 = unique(c(ind1,ind2))
	x = x_in[ind3,]
	}else{
	x = x_in
	ind3 = ind
	}

if(use.loess==TRUE){
	x.loess = loess(read_count~gc,data=x)
	rdcnt_pred = predict(x.loess,newdata=x)
	ind.rm = which(abs(log2((x$read_count+0.001)/rdcnt_pred))>cutoff)

	flag = 0
	k = 1
	max.iter = 100
	while(flag==0&&length(ind.rm)>0&&k<max.iter){
		x.loess = loess(read_count~gc,data=x[-ind.rm,])
		rdcnt_pred1 = predict(x.loess,newdata=x)

		if(max(abs(rdcnt_pred1-rdcnt_pred),na.rm=TRUE)<epsilon){
			flag=1
			}
		rdcnt_pred = rdcnt_pred1
		ind.rm = which(abs(log2((x$read_count+0.001)/rdcnt_pred))>cutoff | is.na(rdcnt_pred))
		k = k+1
		}

	x.fitted= predict(x.loess,newdata=x_in)
	ind.rm = which(abs(log2((x_in$read_count+0.001)/x.fitted))>cutoff | is.na(x.fitted))
	}else{
	x.gam = gam(read_count~s(gc,bs="cr"),data=x)
	ind.rm = which(abs(x.gam$residuals)> 3*sqrt(x.gam$sig2))	
	flag = 0
	k = 1
	max.iter = 100
	x.fitted = x.gam$fitted.values
	while(flag==0&&length(ind.rm)>0&&k<max.iter){
		x.gam = gam(read_count~s(gc,bs="cr"),data=x[-ind.rm,])
		x.fitted1 = predict(x.gam,newdata=x)
		x.residual = x$read_count - x.fitted1
		ind.rm = which(abs(x.residual)> 3*sqrt(x.gam$sig2))
		if(max(x.fitted1-x.fitted,na.rm=TRUE)<epsilon){flag = 1}
		x.fitted = x.fitted1
		k = k+1
		}
	x.fitted = predict(x.gam,newdata=x_in)
	ind.rm = which(abs(x_in$read_count-x.fitted) > 3*sqrt(x.gam$sig2))
	}




ratio = (x_in$read_count+0.001)/x.fitted
if(nrow(x_in)>n.sample){
	ind.sample = sample(ind,n.sample)	
	}else{
	ind.sample = ind
	}
ratio.tmp = ratio[ind.sample]
qt.ratio = quantile(ratio.tmp,probs=c(0.0001,0.9999))
ratio.tmp = ratio.tmp[ratio.tmp>qt.ratio[1] & ratio.tmp<qt.ratio[2]]

ratio.mclust = Mclust(ratio.tmp,modelNames=c("V"))

pvalue.max = sapply





