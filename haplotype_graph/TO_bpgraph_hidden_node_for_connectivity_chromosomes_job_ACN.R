args=commandArgs(T)
t=read.table("SVs.AS_SV.haplotype_phased", stringsAsFactors=F)
ochr=c(1:22, "X");
chr=c();
while(nrow(t)!=0){
	chr=c();
	chr=unique(c(chr,t[1,2],t[1,4]));
	set=c(1);
	i=1;
	while(i <=nrow(t)){
		if(t[i,2] %in% chr){
			if(t[i,4] %in% chr){
				set=unique(c(set,i));
			}else{
				chr=unique(c(chr,t[i,4]));
				set=unique(c(set,i));
				i=1;
				
			}
		}
		if(t[i,4] %in% chr){
			if(t[i,2] %in% chr){
				set=unique(c(set,i));
			}else{
				chr=unique(c(chr,t[i,2]));
				set=unique(c(set,i));
				i=1;
			}	
		}
		i=i+1;
	}
	
	ochr=ochr[!(ochr %in% chr)];

	write.table(t[set,],paste("SVs.AS_SV.haplotype_phased.", paste(sort(chr),collapse="."),sep=""),quote=F, sep="\t",row.names=F, col.names=F);
	system(paste("mkdir ",paste("euler.",paste(sort(chr),collapse="."),sep=""),sep=""));
	system(paste(paste("Rscript ", args[1], "/", "TO_bpgraph_hidden_node_for_connectivity_chromosomes_nohidden_ACN.R ", sep=""), paste(sort(chr),collapse=" ")," ", paste("SVs.AS_SV.haplotype_phased.", paste(sort(chr),collapse="."),sep="")));
	system(paste("cp node_keys degrees edge_information.txt edge_information.txt.nohidden", paste("euler.",paste(sort(chr),collapse="."),sep=""),sep=" "))
	t=t[-set,];
}

if(length(ochr)!=0){
	for( i in 1:length(ochr)){
		chr=ochr[i];
 	       system(paste("mkdir ",paste("euler.",paste(sort(chr),collapse="."),sep=""),sep=""));
 	       system(paste(paste("Rscript ", args[1], "/", "TO_bpgraph_hidden_node_for_connectivity_chromosomes_nohidden_ACN.R ", sep=""), paste(sort(chr),collapse=" ")," ", "NULL",sep=""));
                system(paste("cp node_keys degrees edge_information.txt edge_information.txt.nohidden", paste("euler.",paste(sort(chr),collapse="."),sep=""),sep=" "))
	}
}
