



#infile = "/home/ruibin/software/normalization/normalization/normalization_v0.1.6/tmp/binfiltering_l3MTjXX.txt"

#infile = "/home/ruibin/software/normalization/normalization/normalization_v0.1.6/tmp/binfiltering_7p7rzvv.txt"

#outfile = "/home/ruibin/software/normalization/normalization/normalization_v0.1.6/tmp/binfiltering_copyratio.txt"


args <- commandArgs(TRUE)
if(length(args)!=2)
        { print("Usage: \'R --slave --args <InputData> <copyRatio> < filtering.R\'")
          q()
        }

infile = args[[1]];
outfile = args[[2]];

x_in = read.table(infile,header=T,sep="\t")

ind = c(1:nrow(x_in))
n.sample = 10000
if(nrow(x_in)>n.sample){
	gc.qt = quantile(x_in$gc,probs=c(0.001,0.999))
	ind1 = which(x_in$gc<gc.qt[1] | x_in$gc > gc.qt[2])
	ind2 = sample(ind,n.sample)
	ind3 = unique(c(ind1,ind2)) ## make sure data points with extrem GC contents are always used.
	x = x_in[ind3,]
	}else{
	x = x_in
	ind3 = ind
	}

x.loess = loess(log(read_count+0.1)~gc,data=x)
x.fitted = predict(x.loess,newdata=x_in)
x.fitted = exp(x.fitted)

ratio = (x_in$read_count+0.1)/x.fitted

x.out = data.frame(x_in,ratio=ratio)

write.table(x.out,file=infile,row.names=FALSE,quote=FALSE,sep="\t")

write.table(x.out$ratio,file=outfile,row.names=FALSE,quote=FALSE,sep="\t")


