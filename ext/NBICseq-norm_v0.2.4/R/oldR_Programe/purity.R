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


### now the EM algorithm
## note that mu[3] is mu_2=1 in the manuscript
## length(p) = n+1;
## n is the largest copy number
### E-step function

EM = function(n, x, p0=NULL, alpha0=0.5, tau0=0.2^2, xi0=0,epsilon=1e-6,max.it=100){
	mu = c(0:n)/2
	B = length(x)
	if(is.null(p0)){
		p0 = sapply(1:(n+1), FUN=function(i){return(sum(x<=mu[i]+1/2&x>=mu[i]-1/2)/length(x))})
		}
	
	p_t = p0; alpha_t = alpha0; tau_t = tau0; xi_t = xi0
	p_t1 = p0+1; alpha_t1 = alpha0+1; tau_t1 = tau0+1; xi_t1 = xi0+1
	flag = 0

	rho = matrix(0,nrow=B, ncol=n+1)

	k = 1;
	while(flag==0 && k< max.it){
		### first E-step
		mu_new = mu + alpha_t*(mu[3]-mu) + xi_t
		sqrt.tau_t = sqrt(tau_t)
		for(b in c(1:B)){
			rho.sum = 0
			for(i in c(1:length(mu))){
				rho[b,i] = p_t[i]*dnorm(x[b],mean=mu_new[i],sd=sqrt.tau_t)
				rho.sum = rho.sum + rho[b,i]
				}
			for(i in c(1:length(mu))){
				if(rho.sum!=0) rho[b,i] = rho[b,i]/rho.sum
				}
			}
		
		### then M-step
		a_alpha = 0
		b_alpha = 0
		c_alpha = 0
		c_xi = 0
		for(b in c(1:B)){
			for(i in c(1:length(mu))){
				a_alpha  = a_alpha + rho[b,i]*(mu[3]-mu[i])^2
				b_alpha  = b_alpha + rho[b,i]*(mu[3]-mu[i])
				c_alpha  = c_alpha + rho[b,i]*(mu[3]-mu[i])*(mu[i] - x[b])
				c_xi = c_xi + rho[b,i]*(mu[i] - x[b])
				}
			}
		alpha_t1 = -(B*c_alpha - b_alpha*c_xi)/(B*a_alpha - b_alpha^2)
		xi_t1 = -(a_alpha*c_xi - b_alpha*c_alpha)/(B*a_alpha - b_alpha^2)
		mu_new = mu + alpha_t1*(mu[3]-mu) + xi_t1
		tau_t1 = 0
		for(b in c(1:B)){
			for(i in c(1:length(mu))){
				tau_t1 = tau_t1 + rho[b,i]*(x[b] - mu_new[i])^2
				}
			}
		tau_t1 = tau_t1/B
		p_t1 = colSums(rho)/B

		if(sum(abs(p_t-p_t1))<epsilon && abs(alpha_t1-alpha_t)<epsilon && abs(xi_t1-xi_t)<epsilon && abs(tau_t1-tau_t)<epsilon){flag = 1}
		max.diff = max(sum(abs(p_t-p_t1)),abs(alpha_t1-alpha_t),abs(xi_t1-xi_t),abs(tau_t1-tau_t))
		p_t = p_t1
		alpha_t = alpha_t1
		xi_t = xi_t1
		tau_t = tau_t1
		k = k+1 
		}
	return(list(p=p_t,alpha=alpha_t,xi = xi_t, tau = tau_t))
	}



###
prob = c(0.05,0.1,0.65,0.14,0.05,0.01)
prob = c(0.05,0.1,0.65,0.15,0.05)
mu = c(0:4)/2
sigma2 = rep(0.1^2,length(prob))
alpha = 0.5
xi = 0.1
x = simMixnormal(1000,prob=prob,mu=mu,sigma2=sigma2)
x = (1-alpha)*x + alpha*rnorm(1000,mean=1,sd=sqrt(sigma2[1])) + xi


zz = EM(4,x,alpha=0.1,tau=0.1)
loglk(x,zz$p,zz$alpha,zz$xi,zz$tau) 
loglk(x,prob,alpha,xi,sigma2*((1-alpha)^+alpha^2))

plot(density(x))
xx = c(0:1000)/100
yy = sapply(xx,FUN=function(x){dnormmix(x,zz$p,c(0:(length(zz$p)-1))/2,zz$tau,zz$alpha,zz$xi)})
lines(xx,yy,type="l",xlim=c(0,6),col="red")

yy = sapply(xx,FUN=function(x){dnormmix(x,zz1$p,c(0:(length(zz1$p)-1))/2,zz1$tau,zz1$alpha,zz1$xi)})
lines(xx,yy,type="l",xlim=c(0,6),col="blue")

xx = c(0:1000)/100
yy = sapply(xx,FUN=function(x){dnormmix(x,zz$p,c(0:(length(prob)-1))/2,sigma2[1]*((1-alpha)^2+alpha^2),alpha,xi)})
lines(xx,yy,type="l",xlim=c(0,6),col="blue")




prob = c(0.00,0.0,0.65,0.2,0.15)
mu = c(0:4)
sigma2 = rep(0.1^2,length(prob))
x = simMixnormal(1000,prob=prob,mu=mu,sigma2=sigma2)



#### simulate data from the model as descripted in the manuscript
pi = 0.1
prob = c(0.05,0.1,0.1,0.6,0.1,0.05)
n = length(prob)-1
mu = c(0:n)
D = sum(prob*mu)
sigma2 = 0.06

mu_new = mu*(1-pi) + 2*pi

y = simMixnormal(1000,prob=prob,mu = mu_new, sigma2 = rep(sigma2,length(prob)))

y = y/mean(y)




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




#### try to see the EM results for GCC data
## tumor
dir.copyratio = "/media/Data3/GCC/ploidypurity/ploidypurity/copyratiodata/tumor/"
tumor.files = list.files(dir.copyratio,pattern="^TCGA") 
sample.names = sapply(1:length(tumor.files),FUN=function(i){x=strsplit(tumor.files[i],".",fixed=TRUE); return(x[[1]][1])})
dir.EM = "/media/Data3/GCC/ploidypurity/ploidypurity/EM_estimates/tumor/"
EM.files = paste(sample.names,"_maxComp8.txt",sep="")
plot.dir = "/media/Data3/GCC/ploidypurity/plots/tumor/"

for(s in sample.names){
	copyratio.file = paste(dir.copyratio,s,".tumor.copyratio.txt",sep="")
	EM.file = paste(dir.EM,s,"_maxComp8.txt",sep="")
	plot.file = paste(plot.dir,s,"_tumor_maxComp8.pdf",sep="")
	
	copyratio = scan(copyratio.file)
	em.est = read.purityEM(EM.file)


	y.d = density(copyratio)
	xx = c(0:600)/200
	y.best = sapply(xx,FUN=function(x1){dnormmix(x1,em.est$model.best$p, em.est$model.best$tauv, em.est$model.best$alpha, em.est$model.best$xi)})
	y.bic = sapply(xx,FUN=function(x1){dnormmix(x1,em.est$model.bic$p, em.est$model.bic$tauv, em.est$model.bic$alpha, em.est$model.bic$xi)})

	ylim = c(0,1)
	ylim[2] = max(max(y.d$y),max(y.best),max(y.bic))

	pdf(plot.file,width=7,height=5)
	par(cex=0.9)
	plot(xx,y.best,main=s,col="red",sub="Copy Ratio",ylim=ylim,xlim=c(0,5),type="l",ylab="Density")
	lines(y.d,col="green")
	lines(xx,y.bic,col="blue")
	legend("bottomright",c("Density","Best Model","BIC model"),col=c("green","red","blue"),lty=1)
	l.text = c("Best Model","BIC model")	
	l.text[1] = paste(l.text[1],"\nPloidy: ",em.est$model.best$ploidy,"\nPurity:",em.est$model.best$purity,"\nBIC:",em.est$model.best$bic,"\n")
	l.text[2] = paste(l.text[2],"\nPloidy: ",em.est$model.bic$ploidy,"\nPurity:",em.est$model.bic$purity,"\nBIC:",em.est$model.bic$bic)
	legend("topright",l.text,bty="n")
	dev.off()
	}



	par(cex=0.9,lwd=2)
	plot(xx,y.best,col="red",xlab="Copy Ratio",ylim=ylim,xlim=c(0,5),type="l",ylab="Density")
	lines(y.d,col="blue")
	legend("bottomright",c("Kernal Smoothing","Normal Mixture model"),col=c("blue","red"),lty=1)
	l.text = c("Best Model","BIC model")	
	l.text[1] = paste(l.text[1],"\nPloidy: ",em.est$model.best$ploidy,"\nPurity:",em.est$model.best$purity,"\nBIC:",em.est$model.best$bic,"\n")
	l.text[2] = paste(l.text[2],"\nPloidy: ",em.est$model.bic$ploidy,"\nPurity:",em.est$model.bic$purity,"\nBIC:",em.est$model.bic$bic)
	legend("topright",l.text,bty="n")



