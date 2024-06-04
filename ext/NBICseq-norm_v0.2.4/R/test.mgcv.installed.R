
args <- commandArgs(TRUE)
if(length(args)!=0)
        { print("Usage: \'R --slave < normalize.R\'")
          q()
        }


load.success = library(mgcv,logical.return=TRUE)
if(!load.success){
        q(save="no",status=1)
        }

