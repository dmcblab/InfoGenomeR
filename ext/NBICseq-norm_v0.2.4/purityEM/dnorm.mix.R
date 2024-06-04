

dnorm.mix = function(x,alpha,xi,tau,p){
	n = length(p)
	nu = c(0:(n-1))
	mu = (1-alpha)*nu/2 + alpha + xi
	sdt = sqrt(tau)
	z = sapply(1:n,FUN=function(i){return(p[i]*dnorm(x,mean=mu[i],sd=sdt[i]))})
	return(sum(z))
	}

alpha = 0.322357	
xi = -0.226514	
tau = c(0.00515209,0.0135272,0.0107974,0.00887801,0.0121268,0.00958555)
p = c(0.0390566,0.120674,0.0961385,0.615175,0.0838935,0.0450617)

xx = c(0:1000)/200
y = sapply(1:length(xx),FUN=function(i){return(dnorm.mix(xx[i],alpha,xi,tau,p))})
plot(xx,y,type="l",col="red")

alpha = 0.330772	
xi = -0.221553	
tau = rep(0.00926211,6)
p = c(0.0470017,0.11138,0.0960734,0.620613,0.0770516,0.0478799)


xx = c(0:1000)/200
y = sapply(1:length(xx),FUN=function(i){return(dnorm.mix(xx[i],alpha,xi,tau,p))})
lines(xx,y,type="l",col="blue")



### estimate from test.data.outlier.removed.txt
alpha = 0.242549	
xi = -0.0720961	
tau = c(0.0689529,0.00126059,0.00265714,0.0046158,0.0233873,0.007456,0.00966482)
p = c(0.00383766,0.0373795,0.71497,0.231179,0.00155448,0.0109704,0.000109367)


xx = c(0:1000)/200
y = sapply(1:length(xx),FUN=function(i){return(dnorm.mix(xx[i],alpha,xi,tau,p))})
lines(xx,y,type="l",col="blue")



lk = sapply(1:length(x),FUN=function(i){return(dnorm.mix(x[i],alpha,xi,tau,p))})
log.lk = sum(log(lk))



alpha = 0.122993	
xi = 0.390604	
tau = c(0.0245727,0.0245727,0.0245727,0.0245727,0.0245727)
p = c(0.0125892,0.845986,0.129381,0.0120443,0)


alpha = 0.243188	
xi = -0.0721879	
tau = c(0.0435115,0.00121974,0.00262806,0.0043843,0.249012,32.4472,0,0)
p = c(0.00315749,0.0371798,0.712615,0.226117,0.020522,0.000409104,0,0)




#### estimate from test.data2.txt
alpha = 0.0939244	
xi = -0.0873938	
tau = c(5.0044e-05,0.0850992,0.0177541,0.0493264)
p = c(0.00256484,0.0022286,0.80395,0.191257)


xx = c(0:1000)/200
y = sapply(1:length(xx),FUN=function(i){return(dnorm.mix(xx[i],alpha,xi,tau,p))})
lines(xx,y,type="l",col="blue")

alpha = 0.0928407	
xi = -0.0864868	
tau = c(5.0078e-05,0.0848905,0.0177889,0.0452272,0.278517,0.0714583)
p = c(0.00256375,0.00224464,0.807357,0.187432,1.76147e-06,0.000400703)

xx = c(0:1000)/200
y = sapply(1:length(xx),FUN=function(i){return(dnorm.mix(xx[i],alpha,xi,tau,p))})
lines(xx,y,type="l",col="red")


###TCGA-A6-2671-tumor

alpha = 0.553193	
xi = -0.109339	
tau = c(0.00231515,0.00440502,0.00433637,0.00629152,0.506708,0.0222489)
p = c(0.0132105,0.269696,0.207384,0.36089,0.0264627,0.122357)

xx = c(0:1000)/200
y = sapply(1:length(xx),FUN=function(i){return(dnorm.mix(xx[i],alpha,xi,tau,p))})
lines(xx,y,type="l",col="blue")

alpha = 0.56318	
xi = -0.110065	
tau = c(0.0431759,0.00500426,0.00337906,0.00719299,0.0160286,0.0156117,0.0115826,0.000370807,0.0149372,0.00655684,0.805364,0.659842,0.382499,0.0636999)
p = c(0.023434,0.276589,0.18194,0.385324,0.00176203,0.107606,0.0187034,0.00125879,0.0022984,0.000683459,5.76619e-70,9.98381e-188,4.33625e-72,0.000401307)

xx = c(0:1000)/200
y = sapply(1:length(xx),FUN=function(i){return(dnorm.mix(xx[i],alpha,xi,tau,p))})
lines(xx,y,type="l",col="red")









