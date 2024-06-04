
args <- commandArgs(TRUE)
if(length(args)!=3 && length(args)!=4 )
        { print("Usage: \'R --slave --args <InputData> <output> <fragmentLen> <CopyNumberEstimatFile> < refineGAM.R\'")
          q(save="no",status = 1)
        }


infile = args[[1]];
outfile = args[[2]];
fragLen = as.numeric(args[[3]]);
if(length(args)==4) {cnfile = args[[4]];}

load.success = library(mgcv,logical.return=TRUE)
if(!load.success){
        q(save="no",status=1)
        }


x_in = read.table(infile,header=TRUE);
x_in$ratio = x_in$obs/x_in$expected


x = x_in
ind = c(1:nrow(x))

if(length(args)==4){
	cn = read.table(cnfile,header=TRUE)
	if(nrow(x)!=nrow(cn) || sum(x$start!=cn$start)!=0 || sum(x$end!=cn$end)!=0 || sum(x$obs!=cn$read_count)!=0){
		print("Warning in \'refineGAM.R\': infile and cnfile are not compatible");
		}else{
		cn.table = aggregate(cn$copynumber,by=list(cn$copynumber),length)
		ind.mod = which.max(cn.table[,2])
		ploidy = cn.table[1,ind.mod]
		if(ploidy==2 || cn.table[ind.mod,2]/sum(cn.table[,2])){
			ind = which(cn$copynumber==2)
			}
		}
	}

x = x_in[ind,]
ind1 = which(x$ratio!=0)
z = data.frame(ratio=log(x$ratio+1e-10)[ind1],gc=x$gc[ind1]) ### here ratio is actually log ratio
qnt = quantile(z$ratio, probs=c(.25, .75), na.rm = TRUE) 
H = 2 * IQR(z$ratio, na.rm = TRUE)
z.q = c(qnt[1] - H, qnt[2] + H)
ind2 = which(z$ratio>=z.q[1] & z$ratio<=z.q[2])
z = z[ind2,]

#if(nrow(z)>50000){
#	ind.smpl = sample(1:nrow(z),50000)
#	z = z[ind.smpl,]
#	}


z.model = bam(ratio~s(gc),family=gaussian(),data=z)
z.dta.new = data.frame(ratio=rep(0,2*fragLen+2),gc=c(0:(2*fragLen+1))/(2*fragLen+1))
z.dta.new$ratio = predict(z.model,z.dta.new)


write.table(z.dta.new, file=outfile, row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
