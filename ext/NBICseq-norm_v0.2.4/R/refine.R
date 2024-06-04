
args <- commandArgs(TRUE)
if(length(args)!=3)
        { print("Usage: \'R --slave --args <InputData> <output> <CopyNumberEstimatFile> < refine.R\'")
          q()
        }


infile = args[[1]];
outfile = args[[2]];
cnfile = args[[3]]



x_in = read.table(infile,header=TRUE);
x_in$ratio = x_in$obs/x_in$expected

cn = read.table(cnfile,header=TRUE)

x = x_in
ind = c(1:nrow(x))
if(nrow(x)!=nrow(cn) || sum(x$start!=cn$start)!=0 || sum(x$end!=cn$end)!=0 || sum(x$obs!=cn$read_count)!=0){
	print("Error in \'refine.R\': infile and cnfile are not compatible");
	}else{
	cn.table = aggregate(cn$copynumber,by=list(cn$copynumber),length)
	ind.mod = which.max(cn.table[,2])
	ploidy = cn.table[1,ind.mod]
	if(ploidy==2 || cn.table[ind.mod,2]/sum(cn.table[,2])){
		ind = which(cn$copynumber==2)
		}
	}

x = x_in[ind,]
ind1 = which(x$ratio!=0)
z = data.frame(ratio=log(x$ratio)[ind1],gc=x$gc[ind1])
qnt = quantile(z$ratio, probs=c(.25, .75), na.rm = TRUE) 
H = 2 * IQR(z$ratio, na.rm = TRUE)
z.q = c(qnt[1] - H, qnt[2] + H)
ind2 = which(z$ratio>=z.q[1] & z$ratio<=z.q[2])
z = z[ind2,]


flag = 1
k = 1

z.model = lm(ratio~.,data=z)
coef.new = summary(z.model)$coef
coef = coef.new
while(flag&&k<=2){
	if(max(coef.new[-1,4])>0.05) ## all gc coefficients must be significant
		{flag = 0
		}else{
		coef = coef.new
		k = k+1
		nm = paste("gc",k,sep="")
		z[[nm]]=z$gc^k
		z.model = lm(ratio~.,data=z)
		coef.new = summary(z.model)$coef
		} 
	}

write.table(rbind(coef), file=outfile, row.names=FALSE,col.names=c("estimate","std","zvalue","pvalue"),quote=FALSE,sep="\t")
