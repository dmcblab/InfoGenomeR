args=commandArgs(T);
mode=args[5];
min_size=1000 ### The true set only contains >1000bp somatic SVs.
true=read.table(args[2],stringsAsFactors=F);
true$dif=Inf;
test=read.table(args[1], stringsAsFactors=F);
test=test[test[,2]!=test[,4] | abs(test[,3]-test[,5])>min_size,]
true=true[true[,3] !="X" & true[,6]!="X",]; ## check performance for autochromosomes 
test=test[test[,2]!="X" & test[,4] !="X",];
test$dif=Inf;
break_thres=100;

for(true_i in 1:nrow(true)){
	reverse=strsplit(true[true_i,1],split="to");
	reverse=paste(reverse[[1]][2],"to", reverse[[1]][1],sep="");
	
	m1=which(test$V2==true[true_i,3] & abs(test$V3-true[true_i,4]) < break_thres & test$V4==true[true_i,6] & abs(test$V5-true[true_i,7]) < break_thres & test$V6== true[true_i,1] );
	m2=which(test$V4==true[true_i,3] & abs(test$V5-true[true_i,4]) < break_thres & test$V2==true[true_i,6] & abs(test$V3-true[true_i,7]) < break_thres & test$V6 == reverse)
	
			
	breakdif1=max(abs(true[true_i,4] - test$V3[m1])[1] + 1 , abs(true[true_i,7] - test$V5[m1])[1] + 1);
	breakdif2=max(abs(true[true_i,4] - test$V5[m2])[1] + 1 , abs(true[true_i,7] - test$V3[m2])[1] + 1);


	if(length(m1)!=0){
	       true[true_i,"dif"]=breakdif1;
		test[m1,"dif"]=breakdif1;
	}
	if(length(m2)!=0){
	       true[true_i,"dif"]=breakdif2;
		test[m2,"dif"]=breakdif2;
	}

}


precision=length(which(test$dif<break_thres))/nrow(test)
recall=length(which(true$dif<break_thres))/nrow(true)
Fscore=2*(precision*recall/(precision+recall));

sprintf("precision:%f recall:%f fmeasure: %f",precision, recall, Fscore)
