##### First write a funtion to simulation from mixture of normal models

simMixnormal = function(n, prob, mu=c(1:length(prob)), sigma2=rep(1,length(prob))){
	
	if(length(prob)<1) stop("prob must be a vector of length at least 1")
	if(length(prob)!=length(mu) || length(prob)!=length(sigma2)){stop("prob, mu and sigma2 must be of the same length")}
	if(sum(prob)!=1) stop("Summation of prob must be 1")
	if(sum(prob<0)!=0) stop("prob must be nonnegative numbers")
	if(sum(sigma2<=0)!=0) stop("The variance must be positive")
	sigma = sqrt(sigma2)

	ind = sapply(1:n, FUN=function(i){k = which(as.vector(rmultinom(1,1,prob))!=0); return(k)})

	x = sapply(1:n, FUN=function(i){xx = rnorm(1,mean=mu[ind[i]],sigma[ind[i]]); return(xx)})
	return(x)
	}


dnormmix = function(x,p,mu,tau,alpha,xi){
	if(length(p)!=length(mu)){stop("p and mu must be of the same length")}
	mu_new = mu + alpha*(mu[3]-mu) + xi
	y = sapply(1:length(p),FUN=function(i){p[i]*dnorm(x,mean=mu_new[i],sd=sqrt(tau))})
	return(sum(y))
	}

loglk = function(x,p,alpha,xi,tau){
	mu = c(0:(length(p)-1))/2
	yy = sapply(x,FUN=function(xx){dnormmix(xx,p,mu,tau,alpha,xi)})
	return(sum(log(yy)))
	}



##### functions for reading the purityEM output
read.purityEMone = function(con){
	row1 = read.table(con,header=T,nrows=1)
	tauv = read.table(con,header=T,nrows=1)
	p = read.table(con,header=T,nrows=1)
	if(row1$WF==1){
		psd = read.table(con,header=F,nrows=1)
		pvalue = read.table(con,header=F,nrows=1)
		}else{
		psd = NULL
		pvalue = NULL
		}
	ploidy = as.list(row1)
	ploidy$tauv = as.numeric(tauv)
	ploidy$p = as.numeric(p)
	ploidy$sd.p = as.numeric(psd)
	ploidy$pvalue.p = as.numeric(pvalue)
	return(ploidy)
	}

read.purityEM = function(filename){
	EM = NULL
	con = file(filename,open="r")
	x = readLines(con,n=1)
	i = 1
	while(length(x)>0){
		pushBack(x,con)
		EM[[i]] = read.purityEMone(con)
		i = i+1
		x = readLines(con,n=1)
		}

	bic = sapply(1:length(EM),FUN=function(i){return(EM[[i]]$bic)})
	ploidy = sapply(1:length(EM),FUN=function(i){return(EM[[i]]$ploidy)})
	ind.bic = which.min(bic)
	model.bic = EM[[ind.bic]]
	model.best = EM[[length(EM)]]

	close(con)
	return(list(model.best=model.best,model.bic=model.bic,models.EM=EM,bic=bic,ploidy=ploidy))
	}

dnormmix = function(x,p,tauv,alpha,xi){
        if(length(p)!=length(tauv)){stop("p and tauv must be of the same length")}
	mu = c(1:length(p))/2-1/2
        mu_new = mu + alpha*(mu[3]-mu) + xi
        y = sapply(1:length(p),FUN=function(i){p[i]*dnorm(x,mean=mu_new[i],sd=sqrt(tauv[i]))})
        return(sum(y))
        }

predict.copynumber = function(x,purityModel){
	p = purityModel$p
	tauv = purityModel$tauv
	alpha = purityModel$alpha
	xi = purityModel$xi
	if(length(p)!=length(tauv)){stop("p and tauv must be of the same length")}
	mu = c(1:length(p))/2-1/2
	mu_new = mu + alpha*(mu[3]-mu) + xi
	y = sapply(1:length(p),FUN=function(i){p[i]*dnorm(x,mean=mu_new[i],sd=sqrt(tauv[i]))})
	
	ind = which.max(y)
	#print(y[ind])
	if(y[ind]>1e-6) {
		copynumber = ind - 1
		}else{
		ind = which.min(abs(x-mu_new))	
		copynumber = ind - 1
		if(copynumber == length(mu_new) - 1 ) { copynumber = length(mu_new)}
		}

	return(copynumber)
	}






args <- commandArgs(TRUE)
if(length(args)!=3 && length(args)!=4 && length(args)!=5)
        { print("Usage: \'R --slave --args <copyratioFile> <purityEstimate> <CNVRegionFile> [<FigureFile>] [<FigTitle>] < purity.R \'")
          q()
        }

incopyratio = args[[1]];
inpurityEM = args[[2]];
cnv.region.file = args[[3]];
if(length(args)==4){
	fig.name =args[[4]]
	}
if(length(args)==5){
	fig.title = args[[5]]
	}else{
	fig.title = ""
	}




purity.estimate = read.purityEM(inpurityEM)

x = read.table(incopyratio,header=T)

copynumber.predict = sapply(1:nrow(x),FUN = function(i){predict.copynumber(x$ratio[i],purity.estimate$model.best)})

x$copynumber = copynumber.predict

write.table(x,file=incopyratio,row.names=FALSE,quote=FALSE,sep="\t")

ind = which(copynumber.predict==2) ### nonCNV regions

cnt.mad = mad(x$read_count)
cnt.median = median(x$read_count)
ylim = c(0,cnt.median+10*cnt.mad)


if(length(args)==4 || length(args)==5){
	pdf(fig.name,width=7,height=5)
	par(cex=0.9)
	if(nrow(x)>20000){
		ind.subsmpl = sample(1:nrow(x),20000)
		}else{
		ind.subsmpl = c(1:nrow(x))
		}
	x.subsample = x[ind.subsmpl,]
	plot(x.subsample$gc,x.subsample$read_count,xlab="GC",ylab="Read Count",ylim=ylim,pch=".",col="blue",main=fig.title)
	ind.subsmpl.cnv = which(x.subsample$copynumber!=2)
	if(length(ind.subsmpl.cnv)>0) {points(x.subsample$gc[ind.subsmpl.cnv],x.subsample$read_count[ind.subsmpl.cnv],pch=".",col="red")}

	copyratio = x$ratio
	if(length(copyratio)>10000){
		copyratio = sample(copyratio,10000)
		mdn.cr = median(copyratio)
		mad.cr = mad(copyratio)
		ind = which(copyratio > mdn.cr - 10*mad.cr & copyratio < mdn.cr + 10*mad.cr)
		if(length(ind)>0) copyratio = copyratio[ind]
		}
	y.d = density(copyratio)
	xx = c(0:600)/200
	em.est = purity.estimate
	y.best = sapply(xx,FUN=function(x1){dnormmix(x1,em.est$model.best$p, em.est$model.best$tauv, em.est$model.best$alpha, em.est$model.best$xi)})
	ylim = c(0,1)
	ylim[2] = max(max(y.d$y),max(y.best))
	plot(xx,y.best,col="red",xlab="Copy Ratio",ylim=ylim,xlim=c(0,5),type="l",ylab="Density",main=fig.title)
	lines(y.d,col="blue")
	legend("bottomright",c("Density","Mixture Model"),col=c("blue","red"),lty=1)
	l.text = paste("Ploidy: ",em.est$model.best$ploidy,"\nPurity:",em.est$model.best$purity,"\nBIC:",em.est$model.best$bic,"\n")
	legend("topright",l.text,bty="n")
	dev.off()
	}



if(purity.estimate$model.best$ploidy==2 || purity.estimate$model.best$p[3]>0.2){
	write.table(x[-ind,c(1:3)],file=cnv.region.file,row.names=FALSE,quote=FALSE,sep="\t")
	}else{
	write.table(x[-c(1:nrow(x)),c(1:3)],file=cnv.region.file,row.names=FALSE,quote=FALSE,sep="\t")
	}




