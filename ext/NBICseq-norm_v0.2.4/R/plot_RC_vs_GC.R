
args <- commandArgs(TRUE)
if(length(args)!=2 && length(args)!=3)
        { print("Usage: \'R --slave --args <InputData> <FigFile> < plot_RC_vs_GC.R\'")
	  print("or") 
	  print("Usage: \'R --slave --args <InputData> <FigFile> <Title> < plot_RC_vs_GC.R\'")
          q()
        }



infile = args[[1]];
figfile = args[[2]];

main.title=""
if(length(args)==3){
	main.title = args[[3]]
	}

infiles.name = unlist(strsplit(infile,","))


readcount = c()
gc = c()
for(i in c(1:length(infiles.name))){
	x = read.table(infiles.name[[i]],header=T,sep="\t")
	readcount = c(readcount,x$read_count)
	gc = c(gc,x$gc)
	}


ylim = quantile(readcount,probs=c(0.05,0.95))
yrange = ylim[2] - ylim[1]
ylim[1] = max(ylim[1] - 0.6*yrange,0)
ylim[2] = ylim[2] + 0.6*yrange

pdf(file=figfile,width=7,height=5)
smoothScatter(gc,readcount,xlab="GC",ylab="Read Count",main=main.title,ylim=ylim)
dev.off()
