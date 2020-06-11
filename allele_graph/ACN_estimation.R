library(fitdistrplus)
library(plyr)
args=commandArgs(T)
purity=as.numeric(args[2]);
#NB_MEAN=as.numeric(args[2]);
#NB_DIS=as.numeric(args[3]);
updown=as.numeric(args[5]);
prefiltered_or_not=args[4];
p_thres=0.8
NHEJ_cut_thres=10;
imp_cut_thres=7200;
#args =c();
#args[1]=3;
#library(rtracklayer)
#library(Biostrings)
#library(BSgenome.Hsapiens.UCSC.hg19)

CN=read.table("copy_numbers.CN_opt", header=T, stringsAsFactors=F)
SV=read.table("SVs",stringsAsFactors=F)
SV_opt=read.table("SVs.CN_opt",stringsAsFactors=F)
SV_tag=read.table("SVs.tag",stringsAsFactors=F)
names(SV_tag)[15]="tag";
SV$V2[SV$V2=="X"]=23;
SV$V4[SV$V4=="X"]=23;
SV_tag$V2[SV_tag$V2=="X"]=23;
SV_tag$V4[SV_tag$V4=="X"]=23;
SV_opt$V2[SV_opt$V2=="X"]=23;
SV_opt$V4[SV_opt$V4=="X"]=23;




CN=CN[CN$Chromosome==args[1],];
SV=SV[SV$V2 == args[1] & SV$V4 == args[1],];
SV_tag=SV_tag[SV_tag$V2 == args[1] & SV_tag$V4 == args[1],];

if(args[1] ==23) {
	snp_count=read.table(paste(args[3], "/", "X",".het.count",sep=""), stringsAsFactors=F)
	snp_list=read.table(paste(args[3], "/","X",".het.list",sep=""), stringsAsFactors=F);
}else{
        snp_count=read.table(paste(args[3], "/", args[1],".het.count",sep=""), stringsAsFactors=F)
        snp_list=read.table(paste(args[3], "/", args[1],".het.list",sep=""), stringsAsFactors=F);
}


if(prefiltered_or_not=="T"){
	SNP=read.table("het_snps.format.filtered", stringsAsFactors=F)
}else{
        SNP=read.table("het_snps.format", stringsAsFactors=F)
}


SNP[SNP[,1]=="X",1]<-23

SNP=SNP[SNP$V1==args[1],];

#if(!is.na(SNP[1,7])){
#	SNP[,9]=ceiling(SNP[,9]*SNP[,6]/SNP[,7]);
#}
#SNP=SNP[SNP[,9]!= Inf & !is.nan(SNP[,9]) & !is.na(SNP[,9]),];


#w=which(SNP[,8] > SNP[,9]);
#tmp=SNP[w,8];
#SNP[w,8] = SNP[w,9];
#SNP[w,9] = tmp;


#remove=c();
#for(i in 1:nrow(SNP)){
#	if(length( which(repeatmasker[, (V2 < SNP[i,2] & V3 > SNP[i,2])]))!=0){
#		remove=c(remove,i);
#	}
#}
#if(length(remove)!=0){
#	SNP=SNP[-remove,];
#}

#bigWig="./wgEncodeCrgMapabilityAlign100mer.bigWig";
AS_CN=data.frame();
#AS_CN=data.frame(matrix(,ncol=32,nrow=0),stringsAsFactors=F);
##############
#NB_MEAN=15;
#NB_MEAN_LOWER=as.numeric(NB_MEAN)-updown; ## prevent overfitting
#NB_MEAN_UPPER=as.numeric(NB_MEAN)+updown;
#MAX_BASE_DEPTH=200;
min_segment_length=50;
MAX_BASE_DEPTH=Inf;
MAX_DEPTH_multiple=Inf;
minimum_snp=5;
min_het_snp_density=1e-5;
min_nm_p=0.0001;
min_seg=1e5;
min_conf_seg=1e6;
min_probe=100;
#NB_DIS=30;
MAX_ITER=100;
EM_tol=1e-5;
optim_tol=1e-8;
#purity=0.99;
coeff_IQR=1.5;
#CN_estimate=CN;
#CN_estimate$snp_count=0;
#CN_estimate$mean=0;
#CN_estimate$dis=0;
#for(CN_i in 1:nrow(CN)){
#                t=SNP[SNP$V1==CN[CN_i,1]&SNP$V2>CN[CN_i,2]&SNP$V2<CN[CN_i,3],];
#                t=t[t$V8< MAX_BASE_DEPTH & t$V9< MAX_BASE_DEPTH,];
#                if(nrow(t)> minimum_snp){
#                        t_IQR=IQR(c(t[,8],t[,9]));
#                        Q1=quantile(c(t[,8],t[,9]), 1/4);
#                        Q3=quantile(c(t[,8],t[,9]), 3/4);
#                        t=t[t$V8 > Q1-coeff_IQR*t_IQR & t$V8 < Q3+coeff_IQR*t_IQR & t$V9 > Q1-coeff_IQR*t_IQR & t$V9 < Q3+coeff_IQR*t_IQR,]
#			estimate_hat=fitdist(t$V8+t$V9, "nbinom")$estimate;
#			CN_estimate$snp_count[CN_i] = nrow(t);
#			CN_estimate$mean[CN_i] = estimate_hat[2]/CN[CN_i,18];
#	                CN_estimate$dis[CN_i]=estimate_hat[1];
#		}
#
#}


################
for(CN_i in 1:nrow(CN)){
####################################
	if(CN[CN_i,18]>=1){
#	}else if(CN[CN_i,18]>1){
	#	MAX_DEPTH=CN[CN_i,18]*MAX_DEPTH_multiple*NB_MEAN;
		t=SNP[SNP$V1==CN[CN_i,1]&SNP$V2>CN[CN_i,2]&SNP$V2<CN[CN_i,3],];
		t=t[t$V8< MAX_BASE_DEPTH & t$V9< MAX_BASE_DEPTH,];

		sc_i=min(which(CN[CN_i,2] < snp_count[,2]));
		sc_j=max(which(CN[CN_i,3] > snp_count[,3]));
                p_val=-Inf;
		if(sc_i <= sc_j){
			count=c();
			for(c_i in 4:ncol(snp_count)){
				count[c_i]=sum(as.numeric(snp_count[sc_i:sc_j,c_i]))
			}
			count=sort(count[4:length(count)]);
			count=count[10:(length(count)-10)];
			f=fitdist(count, "nbinom");
			nb_size=f$estimate[1];
			nb_mean=f$estimate[2];
			if(nrow(t)>0){
				lt=t[t[,2] %in% snp_list[,2],];
				p_val=pnbinom(nrow(lt),size=nb_size, mu=nb_mean);
			}else{
                                p_val=pnbinom(0,size=nb_size, mu=nb_mean);
			}
		}

		if(nrow(t)!=0){


	#################### GC MAP FILTER ############
	#
	#		gc_map_filtered_vector=c();
	#		for(t_j in 1:nrow(t)){
	#			######## if t[t_j,2] < 49, error could occur ##########
	#			if(args[1]!=23){
	#				GCfreqM=alphabetFrequency(getSeq(Hsapiens,GRanges(paste("chr",t[t_j,1],sep=""), IRanges(t[t_j,2]-49,t[t_j,2]+49))));
	#			}else {
	#				GCfreqM=alphabetFrequency(getSeq(Hsapiens,GRanges(paste("chr","X",sep=""), IRanges(t[t_j,2]-49,t[t_j,2]+49))));
	#	
	#			}
	#			GCfreq=(GCfreqM[2]+GCfreqM[3])/sum(GCfreqM[1,]);
	#			if(score(import(bigWig, which=GRanges(paste("chr","X",sep=""),t[t_j,2])))>0.6  && GCfreq > 0.3 && GCfreq < 0.6) 
	#				gc_map_filtered_vector=c(gc_map_filtered_vector,t_j);
	#		}
	#		t=t[gc_map_filtered_vector,]
	########################################### 3(IQR) FILTER ###################
			t_IQR=IQR(c(t[,8],t[,9]));
			Q1=quantile(c(t[,8],t[,9]), 1/4);
			Q3=quantile(c(t[,8],t[,9]), 3/4);
			t=t[t$V8 > Q1-coeff_IQR*t_IQR & t$V8 < Q3+coeff_IQR*t_IQR & t$V9 > Q1-coeff_IQR*t_IQR & t$V9 < Q3+coeff_IQR*t_IQR,]
		}

#                t=t[t$V8 + t$V9 < MAX_DEPTH,];

		
############
#                if(nrow(t)>=minimum_snp && nrow(t) / (CN[CN_i,"End.bp"] - CN[CN_i,"Start.bp"]) > min_het_snp_density){
                if(nrow(t)>=minimum_snp && p_val > min_nm_p && nrow(t) < CN[CN_i,"Probes"] && abs(nrow(t) / (CN[CN_i,"End.bp"] - CN[CN_i,"Start.bp"])) > min_het_snp_density){    ####### remove hyper mutation 

#################
                        estimate_hat=fitdist(t$V8+t$V9, "nbinom")$estimate
                        NB_MEAN=estimate_hat[2]/CN[CN_i,18];
                        NB_MEAN_LOWER=max(5,as.numeric(NB_MEAN)-updown); ## prevent overfitting
                        #NB_MEAN_LOWER=as.numeric(NB_MEAN)-updown; ## prevent overfitting
                        NB_MEAN_UPPER=as.numeric(NB_MEAN)+updown;
                        NB_DIS=estimate_hat[1];


			ML_ACN=NA;
			ML_sum=-Inf;
			mixture_L_sum_over_H_sum = 0;
			NB_mean=NB_MEAN;
			NB_dis1=NB_DIS;
			NB_dis2=NB_DIS;


			for(ACN in 0:(floor(CN[CN_i,18]/2))){
				if(ACN!=(CN[CN_i,18]-ACN)){
                #                if(1){

					NB_mean=NB_MEAN;
					NB_dis1=NB_DIS;
					NB_dis2=NB_DIS;
					NB_value=Inf;
					E=matrix(ncol=4, nrow=nrow(t))
					max_iter=MAX_ITER;
					NB_value_diff=Inf;
					while(max_iter ==MAX_ITER || ( max_iter > 0 && NB_value_diff > EM_tol*NB_value)) {
						total_p=dnbinom(t[,8], size=NB_dis1, mu=(ACN*purity+1*(1-purity))*NB_mean)*dnbinom(t[,9], size=NB_dis2, mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*NB_mean)+
							dnbinom(t[,8], size=NB_dis2, mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*NB_mean)*dnbinom(t[,9], size=NB_dis1, mu=(ACN*purity+1*(1-purity))*NB_mean);

						E[,1] = dnbinom(t[,8], size=NB_dis1, mu=(ACN*purity+1*(1-purity))*NB_mean)*dnbinom(t[,9], size=NB_dis2, mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*NB_mean)/total_p;
						E[,2] = 1-E[,1];
						

                                                frr<-function(x){
                                                        sum=0;
							sum=sum(E[,1]*log(dnbinom(t[,8],mu=(ACN*purity+1*(1-purity))*x[1], size=x[2])*dnbinom(t[,9], mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*x[1], size=x[3]))+
								E[,2]*log(dnbinom(t[,8], mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*x[1], size=x[3])*dnbinom(t[,9],mu=(ACN*purity+1*(1-purity))*x[1], size=x[2])))
                                                        -sum;
                                                }
					
						optim_err=0;
						rm(M);

                                        tryCatch({M=optim(c(NB_MEAN,NB_dis1,NB_dis2),frr, method="L-BFGS-B", control=list(factr=optim_tol), lower=c(NB_MEAN_LOWER, 2,2), upper=c(NB_MEAN_UPPER, Inf, Inf))
                                                }, error=function(e){optim_err=1;});

						if(optim_err==0 && exists("M")){
							print(M);
							NB_value_diff=abs(NB_value-M$value);
							NB_mean=M$par[1]
							NB_dis1=M$par[2]
							NB_dis2=M$par[3];
							NB_value=M$value;
##							NB_value = NB_frr(M$par);
							max_iter=max_iter-1;
						}else{
							max_iter=-1;
						}

					}

					mixture_L_sum=0;
					for(m_i in 1:nrow(t)){
						mixture_L_sum = mixture_L_sum + max (dnbinom(t[m_i,8],mu=(ACN*purity+1*(1-purity))*NB_mean, size=NB_dis1)*dnbinom(t[m_i,9], mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*NB_mean, size=NB_dis2),dnbinom(t[m_i,9],mu=(ACN*purity+1*(1-purity))*NB_mean, size=NB_dis1)*dnbinom(t[m_i,8], mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*NB_mean, size=NB_dis2));
					}
							
#					mixture_L_sum=  sum(
#							(max((dnbinom(t[,8],mu=(ACN*purity+1*(1-purity))*NB_mean, size=NB_dis1)*dnbinom(t[,9], mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*NB_mean, size=NB_dis2),
#                                                        dnbinom(t[,9],mu=(ACN*purity+1*(1-purity))*NB_mean, size=NB_dis1)*dnbinom(t[,8], mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*NB_mean, size=NB_dis2))))
#							)

					mixture_L_sum_over_H_sum= mixture_L_sum_over_H_sum + mixture_L_sum;

                                        if( !is.na(mixture_L_sum) &&  ML_sum < mixture_L_sum && NB_mean >=NB_MEAN_LOWER && NB_mean <= NB_MEAN_UPPER){

						ML_ACN=ACN;
						ML_sum=mixture_L_sum;
						ML_mean=NB_mean;
						ML_dis1=NB_dis1;
						ML_dis2=NB_dis2;
                               		 }

				}else{	
					NB_mean=NB_MEAN;
                                        NB_dis=NB_DIS;

                                        frr<-function(x){
                                                sum=0;
						sum=sum(log(dnbinom(t[,8],mu=(ACN*purity+1*(1-purity))*x[1], size=x[2])*dnbinom(t[,9],mu=(ACN*purity+1*(1-purity))*x[1], size=x[2])), na.rm=T);
                                                -sum;
                                       }
					optim_err=0;
                                        rm(M);
                                        tryCatch({M=optim(c(NB_MEAN,NB_dis1,NB_dis2),frr, method="L-BFGS-B", control=list(factr=optim_tol), lower=c(NB_MEAN_LOWER, 2,2), upper=c(NB_MEAN_UPPER, Inf, Inf))
						}, error=function(e){optim_err=1;});
                                        if(optim_err==0 && exists("M")){
						print(M);
						NB_value_diff=abs(NB_value-M$value);
						NB_mean=M$par[1]
						NB_dis=M$par[2]
						NB_value=M$value;
						max_iter=max_iter-1;
					}else{
						max_iter=-1;
					}
                                        mixture_L_sum = 0;
                                        mixture_L_sum=sum((dnbinom(t[,8],mu=(ACN*purity+1*(1-purity))*NB_mean, size=NB_dis)*dnbinom(t[,9],mu=(ACN*purity+1*(1-purity))*NB_mean, size=NB_dis)));

                                        mixture_L_sum_over_H_sum= mixture_L_sum_over_H_sum + mixture_L_sum;

					if( !is.na(mixture_L_sum) &&  ML_sum < mixture_L_sum && NB_mean >=NB_MEAN_LOWER && NB_mean <= NB_MEAN_UPPER){
						ML_ACN=ACN;
						ML_sum=mixture_L_sum;
						ML_mean=NB_mean;
						ML_dis1=NB_dis;
                                                ML_dis2=NB_dis;
					}	

				}
			}

			if(is.na(ML_ACN)!=T){
	                	AS_CN=rbind(AS_CN,cbind(CN[CN_i,],allele_1=ML_ACN,allele_2=CN[CN_i,18]-ML_ACN, balanced=if(ML_ACN==(CN[CN_i,18]-ML_ACN)) 1 else 0, mean=ML_mean, dis1=ML_dis1, dis2=ML_dis2, p=ML_sum/mixture_L_sum_over_H_sum),stringsAsFactors=F);
			}else{
		                AS_CN=rbind(AS_CN,cbind(CN[CN_i,],allele_1=NA,allele_2=NA, balanced=NA, mean=NA, dis1=NA, dis2=NA,p=NA),stringsAsFactors=F);
			}

		}else{
			if(CN[CN_i,3]-CN[CN_i,2] < min_seg){
				AS_CN=rbind(AS_CN,cbind(CN[CN_i,],allele_1=NA,allele_2=NA, balanced=NA, mean=NA, dis1=NA, dis2=NA,p=NA),stringsAsFactors=F);
			}else{
				if( CN[CN_i,3]-CN[CN_i,2] > min_conf_seg){
	                                AS_CN=rbind(AS_CN,cbind(CN[CN_i,],allele_1=CN[CN_i,18],allele_2=0, balanced=0, mean=NA, dis1=NA, dis2=NA,p=1),stringsAsFactors=F);
				}else{
                                        AS_CN=rbind(AS_CN,cbind(CN[CN_i,],allele_1=CN[CN_i,18],allele_2=0, balanced=0, mean=NA, dis1=NA, dis2=NA,p=NA),stringsAsFactors=F);
				}
			}

		}
	


	}else{

        AS_CN=rbind(AS_CN,cbind(CN[CN_i,],allele_1=NA,allele_2=NA, balanced=NA, mean=NA, dis1=NA, dis2=NA,p=NA),stringsAsFactors=F);

	}
		

	print(AS_CN[CN_i,]);
}
AS_CN_backup=AS_CN;

start_z=1;
for( z in 1:nrow(AS_CN)){
        if(start_z &&!is.na( AS_CN[z,"allele_1"]) && AS_CN[z,"allele_1"] != AS_CN[z,"allele_2"]){
                prev_a=AS_CN[z,"allele_1"];
                prev_b=AS_CN[z,"allele_2"];
                start_z=0;
        }else if (!start_z && !is.na(AS_CN[z,"allele_1"]) && AS_CN[z,"allele_1"] != AS_CN[z,"allele_2"]){
		z_diff1_hap = min(1,abs ( prev_a - AS_CN[z,"allele_1"])) + min( 1,abs(prev_b - AS_CN[z,"allele_2"]));
                z_diff1= abs ( prev_a - AS_CN[z,"allele_1"]) + abs(prev_b - AS_CN[z,"allele_2"])
                z_diff2_hap = min(1,abs ( prev_a - AS_CN[z,"allele_2"])) + min( 1,abs(prev_b - AS_CN[z,"allele_1"]));
                z_diff2=  abs ( prev_a - AS_CN[z,"allele_2"]) + abs(prev_b - AS_CN[z,"allele_1"])
                if(z_diff2 < z_diff1 || (z_diff2==z_diff1 && z_diff2_hap < z_diff1_hap)){
                        tmp=AS_CN[z,"allele_1"];
                        AS_CN[z,"allele_1"] = AS_CN[z,"allele_2"];
                        AS_CN[z,"allele_2"] = tmp;
                }
                prev_a=AS_CN[z,"allele_1"];
                prev_b=AS_CN[z,"allele_2"];
        }
}

write.table(AS_CN, paste("copy_numbers.CN_opt.ACN_debug.",args[1],sep=""), sep="\t", quote=F, row.names=F)

if(nrow(AS_CN) == 1 && AS_CN[1,"modal_cn"] == 0){
	AS_CN$ACN_block=0;
	AS_CN[1,"allele_1"]=0;
        AS_CN[1,"allele_2"]=0;
        AS_CN[1,"balanced"]=1;
        write.table(AS_CN, paste("copy_numbers.CN_opt.ACN.",args[1],sep=""), sep="\t", quote=F, row.names=F)
	q();
}


for( i in 1:nrow(AS_CN)){
        if( AS_CN[i,"Probes"] < min_probe){
                AS_CN[i,"allele_1"]=NA;
                AS_CN[i,"allele_2"]=NA;
        }else if(AS_CN[i,"modal_cn"]==0){
                AS_CN[i,"allele_1"]=0;
                AS_CN[i,"allele_2"]=0;
	}
}

pmax=0;
for( i in 1:nrow(AS_CN)){
	if(!is.na(AS_CN[i,"allele_1"]) && !is.na(AS_CN[i,"p"]) && AS_CN[i,"p"] > pmax ){
		pmax=AS_CN[i,"p"];
	}
}
if(pmax > p_thres){
	pmax=p_thres;
}

which_CN_low_conf =  which(is.na(AS_CN$p) | AS_CN$p < pmax);

#if(max(AS_CN$p, na.rm=T) > p_thres){
#	which_CN_low_conf = which(is.na(AS_CN$p) | AS_CN$p < p_thres);
#}else{
#        which_CN_low_conf = which(is.na(AS_CN$p) | AS_CN$p < max(AS_CN$p, na.rm=T));
#}


CN_low_conf = AS_CN[which_CN_low_conf, ];

#CN_low_conf = AS_CN[ AS_CN[,3] - AS_CN[,2] > min_seg & AS_CN[,3] - AS_CN[,2] < min_conf_seg & !is.na(AS_CN[,"allele_1"]), ];
if(nrow(CN_low_conf)>0){
	CN_low_conf$solution_including_this = 0;
}
CN_low_conf <<- CN_low_conf;
#AS_CN[AS_CN[,3] - AS_CN[,2]  < min_conf_seg, "allele_1"]=NA;
#AS_CN[AS_CN[,3] - AS_CN[,2]  < min_conf_seg, "allele_2"]=NA;
#AS_CN[AS_CN[,3] - AS_CN[,2]  < min_conf_seg, "balanced"]=NA;
#AS_CN[AS_CN[,3] - AS_CN[,2]  < min_conf_seg, "mean"]=NA;
#AS_CN[AS_CN[,3] - AS_CN[,2]  < min_conf_seg, "dis1"]=NA;
#AS_CN[AS_CN[,3] - AS_CN[,2]  < min_conf_seg, "dis2"]=NA;
AS_CN[which_CN_low_conf, "allele_1"]=NA;
AS_CN[which_CN_low_conf, "allele_2"]=NA;
AS_CN[which_CN_low_conf, "balanced"]=NA;
AS_CN[which_CN_low_conf, "mean"]=NA;
AS_CN[which_CN_low_conf, "dis1"]=NA;
AS_CN[which_CN_low_conf, "dis2"]=NA;





AS_CN_opt = AS_CN;
min_opt=Inf;
min_hap1_diff=Inf;
min_hap2_diff=Inf;
min_hap_diff = Inf;
solution_count=0;

recursive_opt <- function( r_start, r_end, r_i, l_AS_CN, is_greedy){
	
	if(is_greedy && r_i <= r_end){
		a1=min ( l_AS_CN[r_start-1,"allele_1"], l_AS_CN[r_start-1,"allele_2"], l_AS_CN[r_end+1,"allele_1"], l_AS_CN[r_end+1,"allele_2"], na.rm=T);

		l_AS_CN[r_i,"allele_1"] = a1;
		l_AS_CN[r_i,"allele_2"] = max(0, l_AS_CN[r_i,"modal_cn"] - l_AS_CN[r_i,"allele_1"]);
        #	print(r_i)
	        if(l_AS_CN[r_i,"allele_1"] ==  l_AS_CN[r_i,"allele_2"]){
                        l_AS_CN[r_i,"balanced"]=1;
                }else{
                        l_AS_CN[r_i,"balanced"]=0;
                }
                recursive_opt(r_start,r_end,r_i+1, l_AS_CN, is_greedy);

	}else{
		for( v in 0:l_AS_CN[r_i,"modal_cn"]){
			#print("recur");
			l_AS_CN[r_i,"allele_1"]=v;
			l_AS_CN[r_i,"allele_2"]=l_AS_CN[r_i,"modal_cn"]-v;
			if( r_i < r_end){
					recursive_opt(r_start, r_end, r_i+1, l_AS_CN, is_greedy);
			}else{
				iter=0;
				while(iter<=1){
					if(iter==1){
						if(r_start!=1){
							tmp=l_AS_CN[r_start-1,"allele_1"];
							l_AS_CN[r_start-1,"allele_1"]=l_AS_CN[r_start-1,"allele_2"]
							l_AS_CN[r_start-1,"allele_2"]=tmp;
						}else{

							tmp=l_AS_CN[r_end+1,"allele_1"];
							l_AS_CN[r_end+1,"allele_1"]=l_AS_CN[r_end+1,"allele_2"]
							l_AS_CN[r_end+1,"allele_2"]=tmp;
						}
					}
			

					diff_val=0;
					hap1_diff=0;
					hap2_diff=0;
					for(l_AS_CN_li in r_start:r_end){
						if(l_AS_CN_li != 1){
							hap1_diff = hap1_diff + min(1, abs(l_AS_CN[l_AS_CN_li-1,"allele_1"] - l_AS_CN[l_AS_CN_li, "allele_1"]));
							hap2_diff = hap2_diff + min(1, abs(l_AS_CN[l_AS_CN_li-1,"allele_2"] - l_AS_CN[l_AS_CN_li, "allele_2"]));
							diff_val= diff_val + abs(l_AS_CN[l_AS_CN_li-1,"allele_1"] - l_AS_CN[l_AS_CN_li, "allele_1"]) + abs(l_AS_CN[l_AS_CN_li-1,"allele_2"] - l_AS_CN[l_AS_CN_li, "allele_2"])
						}
					}
					if(l_AS_CN_li !=nrow(AS_CN)){
						hap1_diff = hap1_diff + min(1, abs(l_AS_CN[l_AS_CN_li,"allele_1"] - l_AS_CN[l_AS_CN_li+1, "allele_1"]));
						hap2_diff = hap2_diff + min(1, abs(l_AS_CN[l_AS_CN_li,"allele_2"] - l_AS_CN[l_AS_CN_li+1, "allele_2"]));
						diff_val = diff_val + abs(l_AS_CN[l_AS_CN_li,"allele_1"] - l_AS_CN[l_AS_CN_li+1, "allele_1"]) + abs(l_AS_CN[l_AS_CN_li,"allele_2"] - l_AS_CN[l_AS_CN_li+1, "allele_2"])
					}

		
					hap_diff = min(1,hap1_diff) + min(1, hap2_diff);				

 #                                       if(diff_val < min_opt ||
  #                                        (diff_val == min_opt && hap1_diff + hap2_diff < min_hap1_diff + min_hap2_diff ) ||
   #                                       (diff_val == min_opt && hap1_diff + hap2_diff == min_hap1_diff + min_hap2_diff && hap_diff <= min_hap_diff )){
#						if(diff_val == min_opt && hap1_diff + hap2_diff == min_hap1_diff + min_hap2_diff && hap_diff  == min_hap_diff){

                                        if(diff_val < min_opt ||
                                          (diff_val == min_opt &&  hap_diff <= min_hap_diff )){
                                                if(diff_val == min_opt && hap_diff  == min_hap_diff){

#							print( l_AS_CN[l_AS_CN_li,"allele_1"]);
 #                                                       print( l_AS_CN[l_AS_CN_li,"allele_2"]);
#							print(diff_val);
							if(!((AS_CN_opt[l_AS_CN_li,"allele_1"] == l_AS_CN[l_AS_CN_li,"allele_1"] && AS_CN_opt[l_AS_CN_li,"allele_2"] == l_AS_CN[l_AS_CN_li,"allele_2"]) ||
							(AS_CN_opt[l_AS_CN_li,"allele_1"] == l_AS_CN[l_AS_CN_li,"allele_2"] && AS_CN_opt[l_AS_CN_li,"allele_2"] == l_AS_CN[l_AS_CN_li,"allele_1"]))){
								solution_count<<-solution_count + 1;
							}
						}else{
#							print(diff_val);
							solution_count<<-1;
						}
						min_opt<<-diff_val;
						min_hap1_diff <<- hap1_diff;
						min_hap2_diff <<- hap2_diff;
						min_hap_diff <<- hap_diff;
						AS_CN_opt <<- l_AS_CN;
					}

					iter=iter+1;
				}
			}
		}
	}
}

if(nrow(AS_CN)!=1){
	i=1;
	while(i <=nrow(AS_CN)){
                if(is.na(AS_CN[i,"allele_1"] )){
	#	if(is.na(AS_CN[i,"allele_1"] ) && AS_CN[i,"modal_cn"]!=0){
			j=i;
                        while(j+1 <= nrow(AS_CN) && is.na(AS_CN[j+1, "allele_1"])){
			#while(j+1 <= nrow(AS_CN) && is.na(AS_CN[j+1, "allele_1"]) && AS_CN[j+1, "modal_cn"]!=0){
				j=j+1;
			}
			min_opt=Inf;
			min_hap1_diff=Inf;
			min_hap2_diff=Inf;
			min_hap_diff = Inf;
			solution_count=0;

	##		if(prod(AS_CN$modal_cn[i:j])> 1e3){
			if(i!=j){
			#	recursive_opt(i,j,i, AS_CN,T);
			}else{
                  ##              recursive_opt(i,j,i, AS_CN,F);
			}
			
			if(solution_count == 1){
				AS_CN=AS_CN_opt;
				for(AS_CN_i in i:j){
					if(AS_CN[AS_CN_i,"allele_1"]==AS_CN[AS_CN_i,"allele_2"]){
						AS_CN[AS_CN_i,"balanced"]=1;
					}else{
						AS_CN[AS_CN_i,"balanced"]=0;
					}
				}
			}
			i=j+1;
		}else{
			i=i+1;
		}
	}
}



imp = data.frame(matrix(,nrow=0,ncol=6));
imp_i=1;
colnames(imp)=c("start","end","c1","c2","c3","c4");
if(nrow(AS_CN)>  1){
	i=1;
	while(i <= nrow(AS_CN)){
                if(is.na(AS_CN[i,"allele_1"] ) ){
                        j=i;
                        while(j+1 <= nrow(AS_CN) && is.na(AS_CN[j+1, "allele_1"])){
                                j=j+1;
                        }
			c_candidate=c(AS_CN[i-1,"allele_1"],AS_CN[i-1,"allele_2"], AS_CN[j+1, "allele_1"], AS_CN[j+1, "allele_2"]);
			c_candidate=c_candidate[!is.na(c_candidate)];

			count_v=c(AS_CN_backup[i:j, "allele_1"] , AS_CN_backup[i:j, "allele_2"]);
			count_v= count_v[!is.na(count_v)];
			if(length(count_v) > 0){
				count_a= count(count_v);
				count_a=count_a[order(count_a$freq),];
			}

#                        if(length(unique(c_candidate[duplicated(c_candidate)])) == 1 ){
			if(length(unique(c_candidate[duplicated(c_candidate)])) == 1 && ((length(count_v)==0) || ( length(count_v) > 0 && count_a[nrow(count_a),1] == unique(c_candidate[duplicated(c_candidate)]))) ){
				c_candidate = unique(c_candidate[duplicated(c_candidate)]);
			}else{
				c_candidate=unique(c_candidate[!is.na(c_candidate)]);
			}
			imp[imp_i,"start"]=i;
                        imp[imp_i,"end"]=j;
			if(length(c_candidate)>0){
				for(c_i in 1:length(c_candidate)){
					imp[imp_i,2+c_i]=c_candidate[c_i];
				}
			}
			imp_i=imp_i+1;
			i=j+1;	
                }else{
			i=i+1;
		}
        }
}


recursive_NHEJ_opt_AS_CN=AS_CN;

SVs_max_sum=0;

recursive_NHEJ_opt <- function(l_AS_CN,l_opts_unit, l_SV,r_i, l_SV_idx1, l_SV_idx2, l_SV_idx3, l_SV_idx4){
	if( r_i <= length(l_opts_unit)){
		recursive_NHEJ_opt(l_AS_CN, l_opts_unit, l_SV, r_i+1, l_SV_idx1, l_SV_idx2, l_SV_idx3, l_SV_idx4);
		tmp=l_AS_CN[which(l_AS_CN$ACN_block == l_opts_unit[r_i]),"allele_1"];
		l_AS_CN[which(l_AS_CN$ACN_block == l_opts_unit[r_i]),"allele_1"]=l_AS_CN[which(l_AS_CN$ACN_block == l_opts_unit[r_i]),"allele_2"];
		l_AS_CN[which(l_AS_CN$ACN_block == l_opts_unit[r_i]),"allele_2"]=tmp;
		recursive_NHEJ_opt(l_AS_CN, l_opts_unit, l_SV, r_i+1, l_SV_idx1, l_SV_idx2, l_SV_idx3, l_SV_idx4);
	}else{

		SVs_sum=0;

		v1=abs (l_AS_CN[l_SV_idx1,"allele_1"] - l_AS_CN[l_SV_idx2,"allele_1"])
		v2=abs (l_AS_CN[l_SV_idx3,"allele_1"] - l_AS_CN[l_SV_idx4,"allele_1"])
		w1=which(v1>v2);
		if(length(w1) >0){
			v1[w1] = v2[w1];
		}

		v3=abs ( l_AS_CN[l_SV_idx1,"allele_2"] - l_AS_CN[l_SV_idx2,"allele_2"])
		v4=abs ( l_AS_CN[l_SV_idx3,"allele_2"] - l_AS_CN[l_SV_idx4,"allele_2"])
		w2=which(v3>v4);
		if(length(w2) >0){
			v3[w2] = v4[w2];
		}
		w3=which(v1<v3);
		if(length(w3) > 0){
			v1[w3] = v3[w3];
		}
		SVs_sum = sum(v1, na.rm=T);


		if(SVs_max_sum <= SVs_sum){
		       SVs_max_sum <<- SVs_sum;
			recursive_NHEJ_opt_AS_CN <<- l_AS_CN;
		}
	}
}


return_block = function ( l_SV, l_AS_CN){
	ori=strsplit(l_SV[1,6], split="to")[[1]];
	if(ori[1] == "3"){
		w1=which(l_AS_CN$End.bp == l_SV[1,3]);
		w2=which(l_AS_CN$Start.bp == l_SV[1,3]+1);
		block1=unique(l_AS_CN[l_AS_CN$End.bp == l_SV[1,3],"ACN_block"], l_AS_CN[l_AS_CN$Start.bp == l_SV[1,3]+1,"ACN_block"])[1];
	}else{
                w1=which(l_AS_CN$Start.bp == l_SV[1,3]);
                w2=which(l_AS_CN$End.bp == l_SV[1,3]-1);
                block1=unique(l_AS_CN[l_AS_CN$Start.bp == l_SV[1,3],"ACN_block"], l_AS_CN[l_AS_CN$End.bp == l_SV[1,3]-1,"ACN_block"])[1];
	}
        if(ori[2] == "3"){
                w3=which(l_AS_CN$End.bp == l_SV[1,5]);
                w4=which(l_AS_CN$Start.bp == l_SV[1,5]+1);
                block2=unique(l_AS_CN[l_AS_CN$End.bp == l_SV[1,5],"ACN_block"], l_AS_CN[l_AS_CN$Start.bp == l_SV[1,5]+1,"ACN_block"])[1];
        }else{
                w3=which(l_AS_CN$Start.bp == l_SV[1,5]);
                w4=which(l_AS_CN$End.bp == l_SV[1,5]-1);
                block2=unique(l_AS_CN[l_AS_CN$Start.bp == l_SV[1,5],"ACN_block"], l_AS_CN[l_AS_CN$End.bp == l_SV[1,5]-1,"ACN_block"])[1];
        }
	return (c(block1,block2,w1,w2,w3,w4))
}


r_diff=Inf;
min_edit_d_from_low_conf=Inf;
SV_sum_max=0;
g_length=0;
AS_CN_opt=AS_CN;
solution_count=0;
recursive_imp = function(l_AS_CN,r_i, r_end, l_imp){
	if(r_i <= r_end){
		for(i in 3:ncol(l_imp)){	
			if(!is.na(l_imp[r_i,i])){
				w=l_imp[r_i,"start"]:l_imp[r_i,"end"];
				l_AS_CN[w,"allele_1"] = l_imp [r_i,i];
				w1=which(l_AS_CN$modal_cn<l_AS_CN$allele1 & l_AS_CN[,3] - l_AS_CN[,2] > min_segment_length);
				if(length(w1) >0){
					l_AS_CN[w1,"allele_1"] =l_AS_CN[w1,"modal_cn"];
				}

				l_AS_CN[w,"allele_2"] = l_AS_CN[w,"modal_cn"] - l_AS_CN[w,"allele_1"];
				w2= which(l_AS_CN[,"allele_2"]  <0);
				if(length(w2) >0){
					l_AS_CN[w2, "allele_2"] = 0;
				}

				recursive_imp (l_AS_CN, r_i+1, r_end, l_imp);
			}
		}
	}else{
		l_AS_CN[is.na(l_AS_CN[,"allele_1"]),"balanced"] = NA;
		l_AS_CN[which(l_AS_CN$allele_1 == l_AS_CN$allele_2),"balanced"] = 1;
		l_AS_CN[which(l_AS_CN$allele_1 != l_AS_CN$allele_2),"balanced"] = 0;

	
		if(nrow(l_AS_CN) > 1 && length(which(!is.na(l_AS_CN[,"allele_1"]))) > 0  ){
			for(i in 1:(nrow(l_AS_CN)-1)){
				
				SV_flanked=0;
				 w=which  ( (SV_opt[,2] == l_AS_CN[i,"Chromosome"]  & SV_opt[,3] == l_AS_CN[i,"End.bp"] ) | (SV_opt[,4] == l_AS_CN[i,"Chromosome"]  & SV_opt[,5] == l_AS_CN[i,"End.bp"] )  |
					   (SV_opt[,2] == l_AS_CN[i+1,"Chromosome"]  & SV_opt[,3] == l_AS_CN[i+1,"Start.bp"] ) | (SV_opt[,4] == l_AS_CN[i+1,"Chromosome"]  & SV_opt[,5] == l_AS_CN[i+1,"Start.bp"] )
					);

				if(length(w) !=0){
					SV_flanked=1;
					ini_SV_CN=SV_opt[w,15][1];
				}

	#			if(!is.na(l_AS_CN[i,"balanced"]) && !is.na(l_AS_CN[i+1,"balanced"]) && l_AS_CN[i,"balanced"]==0 && l_AS_CN[i+1,"balanced"]==0 && l_AS_CN[i,1]==l_AS_CN[i+1,1]){
                                if(!is.na(l_AS_CN[i,"balanced"]) && !is.na(l_AS_CN[i+1,"balanced"]) && l_AS_CN[i,1]==l_AS_CN[i+1,1]){

					ACN_hap_diff1=min(1,abs(l_AS_CN[i,"allele_1"]-l_AS_CN[i+1,"allele_1"]))  + min(1, abs(l_AS_CN[i,"allele_2"]-l_AS_CN[i+1,"allele_2"]));
					ACN_diff1=abs(l_AS_CN[i,"allele_1"]-l_AS_CN[i+1,"allele_1"])+abs(l_AS_CN[i,"allele_2"]-l_AS_CN[i+1,"allele_2"]);
					ACN_min1=min(abs(l_AS_CN[i,"allele_1"]-l_AS_CN[i+1,"allele_1"]),abs(l_AS_CN[i,"allele_2"]-l_AS_CN[i+1,"allele_2"]));

					ACN_hap_diff2=min(1,abs(l_AS_CN[i,"allele_1"]-l_AS_CN[i+1,"allele_2"]))+min(1,abs(l_AS_CN[i,"allele_2"]-l_AS_CN[i+1,"allele_1"]));
					ACN_diff2=abs(l_AS_CN[i,"allele_1"]-l_AS_CN[i+1,"allele_2"])+abs(l_AS_CN[i,"allele_2"]-l_AS_CN[i+1,"allele_1"]);
					ACN_min2=min(abs(l_AS_CN[i,"allele_1"]-l_AS_CN[i+1,"allele_2"]),abs(l_AS_CN[i,"allele_2"]-l_AS_CN[i+1,"allele_1"]));


                                        if(SV_flanked){
                                                ACN_diff1=abs(ACN_diff1-ini_SV_CN);
						ACN_diff2=abs(ACN_diff2-ini_SV_CN);

					}

					if((ACN_diff2 < ACN_diff1) ||
					   (ACN_diff2 == ACN_diff1 && ACN_min2 < ACN_min1) ||
					   (ACN_diff2 == ACN_diff1 && ACN_min2 == ACN_min1 &&  ACN_hap_diff2 <=ACN_hap_diff1)){
						temp=l_AS_CN[i+1,"allele_1"];
						l_AS_CN[i+1,"allele_1"]=l_AS_CN[i+1,"allele_2"];
						l_AS_CN[i+1,"allele_2"]=temp;
					}

				}
			}
		}
 
		l_AS_CN$ACN_block=NA;
		ACN_block=0;
		block_start=0;
		for( i in 1:nrow(l_AS_CN)){
#			if(l_AS_CN[i,"modal_cn"]==0){
#				l_AS_CN[i,"allele_1"]=0;
#				l_AS_CN[i,"allele_2"]=0;
#			}

			if(!is.na(l_AS_CN[i,"balanced"]) && l_AS_CN[i,"balanced"] == 0 ){
				l_AS_CN[i,"ACN_block"] = ACN_block;
				block_start=1;
			}else{
				if(block_start==1){
					ACN_block=ACN_block+1;
					block_start=0;
				}
			}
		}


		new_ACN_block=0;
		NHEJ_opts= unique(l_AS_CN$ACN_block[complete.cases(l_AS_CN$ACN_block)]);
		test_l_AS_CN <<- l_AS_CN;
		while(length(NHEJ_opts)!=0){
			opts_unit = NHEJ_opts[1];
			
			if(length(which(l_AS_CN$ACN_block == opts_unit)) != 0){
				p_opts_unit=c();
				while(length(p_opts_unit) != length(opts_unit)){
					p_opts_unit = opts_unit;
					nodes = c();

					test<<-l_AS_CN$ACN_block;
					test_opts_unit<<-opts_unit;


					for( j in 1:length(opts_unit)){
						where_block = which(l_AS_CN$ACN_block == opts_unit[j]);
						nodes = c(nodes, l_AS_CN[where_block, 2],  l_AS_CN[where_block, 3]);
						if(where_block[1] != 1){
							nodes= c(nodes, l_AS_CN[where_block[1] -1 ,3]);
						}
						if(where_block[length(where_block)] != nrow(l_AS_CN)){
							nodes = c(nodes, l_AS_CN[where_block[length(where_block)]+1, 2]);
						}
					}
					l_SV = SV[(SV$V3 %in% nodes | SV$V5 %in% nodes ) & SV_tag$tag == 0 , ] ;
#					SV_nodes = c( SV[SV$V3 %in% nodes, "V3"],  SV[SV$V3 %in% nodes, "V5"], SV[SV$V5 %in% nodes, "V3"], SV[SV$V5 %in% nodes, "V5"]);
					SV_nodes = c(l_SV$V3, l_SV$V5);
					opts_unit = c(opts_unit, l_AS_CN[l_AS_CN[,2] %in% SV_nodes,"ACN_block"], l_AS_CN[l_AS_CN[,3] %in% (SV_nodes-1), "ACN_block"], l_AS_CN[l_AS_CN[,3] %in% SV_nodes,"ACN_block"] , l_AS_CN[l_AS_CN[,2] %in% (SV_nodes+1), "ACN_block"])

					opts_unit = unique(opts_unit[complete.cases(opts_unit)]);
				}
				SVs_max_sum<<-0;
				if(nrow(l_SV)>0){
                                        recursive_NHEJ_opt_AS_CN<<-l_AS_CN;
					l_SV_idx1=c();
					l_SV_idx2=c();
					l_SV_idx3=c();
					l_SV_idx4=c();
					for(r_j in 1:nrow(l_SV)){
						ori=strsplit(l_SV[r_j,6], split="to")[[1]];
						if(ori[1] == "5"){
							l_SV_idx1=c(l_SV_idx1,which(l_AS_CN[,2] == l_SV[r_j,3]));
							l_SV_idx2=c(l_SV_idx2,which(l_AS_CN[,3] == l_SV[r_j,3]-1));
						}else{
							l_SV_idx1 = c(l_SV_idx1, which(l_AS_CN[,3] == l_SV[r_j,3]));
							l_SV_idx2 = c(l_SV_idx2, which(l_AS_CN[,2] == l_SV[r_j,3]+1));
						}
						if(ori[2] == "5" ){
                                                        l_SV_idx3=c(l_SV_idx3,which(l_AS_CN[,2] == l_SV[r_j,5]));
                                                        l_SV_idx4=c(l_SV_idx4,which(l_AS_CN[,3] == l_SV[r_j,5]-1));
						}else{
                                                        l_SV_idx3 = c(l_SV_idx3, which(l_AS_CN[,3] == l_SV[r_j,5]));
                                                        l_SV_idx4 = c(l_SV_idx4, which(l_AS_CN[,2] == l_SV[r_j,5]+1));
						}
					}

                                        opts_unit_solve= opts_unit;
                                        while(length(opts_unit_solve) > 0){
                                                opts_unit_cut = opts_unit_solve[1: min(NHEJ_cut_thres, length(opts_unit_solve))];
                                                recursive_NHEJ_opt(l_AS_CN, opts_unit_cut, l_SV, 1, l_SV_idx1, l_SV_idx2, l_SV_idx3, l_SV_idx4);
                                                l_AS_CN=recursive_NHEJ_opt_AS_CN;
                                                opts_unit_solve = opts_unit_solve [-(1:length(opts_unit_cut))];
                                        }
				}

				
				l_AS_CN[l_AS_CN$ACN_block %in% opts_unit, "ACN_block"] = new_ACN_block;
				new_ACN_block=new_ACN_block+1;
			}

			NHEJ_opts = NHEJ_opts[!(NHEJ_opts %in% opts_unit)];
		}




		t=l_AS_CN;
		d=list();
		for( d_i in 1:2){
			d[[d_i]]=data.frame(chrom=t[,1],key=0, degree=0, stringsAsFactors=F)
			dindex=1
			for (i in 1:nrow(t)){
				d[[d_i]][dindex,1]=t[i,1]
				d[[d_i]][dindex,2]=t[i,2]
				d[[d_i]][dindex,3]=t[i,paste("allele_",d_i,sep="")]
				dindex=dindex+1;
				d[[d_i]][dindex,1]=t[i,1]
				d[[d_i]][dindex,2]=t[i,3]
				d[[d_i]][dindex,3]=t[i,paste("allele_",d_i,sep="")]
				dindex=dindex+1;
			 }
		}
		
		if(nrow(SV_opt)>0){
			for( i in 1:nrow(SV_opt)){
				if(SV_opt[i,2] != SV_opt[i,4]){
					if(SV_opt[i,2] == l_AS_CN[1,"Chromosome"]){
						ori=strsplit(SV_opt[i,6],split="to")[[1]];
						if(ori[1] == "5"){
							av_d=-1;
						}else{
							av_d=1;
						}
						av1=abs (d[[1]][ d[[1]][,2] == SV_opt[i,3], 3] - SV_opt[i,15] - d[[1]][ which(d[[1]][,2] == SV_opt[i,3]) + av_d , 3] );
						av2=abs (d[[2]][ d[[2]][,2] == SV_opt[i,3], 3] - SV_opt[i,15] - d[[2]][ which(d[[2]][,2] == SV_opt[i,3]) + av_d , 3] );
						if(!is.na(av1) && !is.na(av2)){
							if(av1 < av2){
								d[[1]][ d[[1]][,2] == SV_opt[i,3], 3] = d[[1]][ d[[1]][,2] == SV_opt[i,3], 3] - SV_opt[i,15];
							}else{
								d[[2]][ d[[2]][,2] == SV_opt[i,3], 3] = d[[2]][ d[[2]][,2] == SV_opt[i,3], 3] - SV_opt[i,15];
							}
						}
					}
                                        if(SV_opt[i,4] == l_AS_CN[1,"Chromosome"]){
                                                ori=strsplit(SV_opt[i,6],split="to")[[1]];
                                                if(ori[2] == "5"){
                                                        av_d=-1;
                                                }else{
                                                        av_d=1;
                                                }
                                                av1=abs (d[[1]][ d[[1]][,2] == SV_opt[i,5], 3] - SV_opt[i,15] - d[[1]][ which(d[[1]][,2] == SV_opt[i,5]) + av_d , 3] );
                                                av2=abs (d[[2]][ d[[2]][,2] == SV_opt[i,5], 3] - SV_opt[i,15] - d[[2]][ which(d[[2]][,2] == SV_opt[i,5]) + av_d , 3] );
						if(!is.na(av1) && !is.na(av2)){
	                                                if(av1 < av2){
	                                                        d[[1]][ d[[1]][,2] == SV_opt[i,5], 3] = d[[1]][ d[[1]][,2] == SV_opt[i,5], 3] - SV_opt[i,15];
	                                                }else{
	                                                        d[[2]][ d[[2]][,2] == SV_opt[i,5], 3] = d[[2]][ d[[2]][,2] == SV_opt[i,5], 3] - SV_opt[i,15];
	                                                }
						}
                                        }
				}
			}
		}

#		l_AS_CN_degree = l_AS_CN;
		SV_sum=0;
		if(nrow(SV)>0){
			for( i in 1:nrow(SV)){
				ori=strsplit(SV[i,6],split="to")[[1]];
				if(ori[1] == 3){
					i1=which(l_AS_CN[,"End.bp"]==SV[i,3]);
					i2=3;
				}else{
					i1=which(l_AS_CN[,"Start.bp"]==SV[i,3]);
					i2=2;
				}
				if(ori[2] == 3){
					i3=which(l_AS_CN[,"End.bp"]==SV[i,5]);
					i4=3;
				}else{
					i3=which(l_AS_CN[,"Start.bp"]==SV[i,5]);
					i4=2;
				}
                                if( !is.na(l_AS_CN[i1, "ACN_block"]) && !is.na(l_AS_CN[i3, "ACN_block"]) && l_AS_CN[i1, "ACN_block"]  != l_AS_CN[i3, "ACN_block"] ){
##				if( is.na(l_AS_CN[i1, "ACN_block"]) || is.na(l_AS_CN[i3, "ACN_block"]) || l_AS_CN[i1, "ACN_block"]  != l_AS_CN[i3, "ACN_block"] ){
					fi1=which.max(c(max(0, l_AS_CN[i1, "allele_1"] - l_AS_CN[i1 + if(i2==3) 1 else -1  , "allele_1"]) , max(0, l_AS_CN[i1, "allele_2"] - l_AS_CN[i1 + if(i2==3) 1 else -1  , "allele_2"]) ));
					fi2=which.max(c(max(0, l_AS_CN[i3, "allele_1"] - l_AS_CN[i3 + if(i4==3) 1 else -1  , "allele_1"]) , max(0, l_AS_CN[i3, "allele_2"] - l_AS_CN[i3 + if(i4==3) 1 else -1  , "allele_2"]) ));

					sv_cn=min ( max (max(0, l_AS_CN[i1, "allele_1"] - l_AS_CN[i1 + if(i2==3) 1 else -1  , "allele_1"]) , max(0, l_AS_CN[i1, "allele_2"] - l_AS_CN[i1 + if(i2==3) 1 else -1  , "allele_2"]) ) ,
						    max (max(0, l_AS_CN[i3, "allele_1"] - l_AS_CN[i3 + if(i4==3) 1 else -1  , "allele_1"]) , max(0, l_AS_CN[i3, "allele_2"] - l_AS_CN[i3 + if(i4==3) 1 else -1  , "allele_2"])) );
				
					if(!is.na(sv_cn)){
						d[[fi1]][ 2*(i1-1)+1 +  if(i2==3) 1 else 0 , "degree"] = d[[fi1]][ 2*(i1-1)+1 +  if(i2==3) 1 else 0 , "degree"] - sv_cn;
	                                       d[[fi2]][ 2*(i3-1)+1 +  if(i4==3) 1 else 0 , "degree"] = d[[fi2]][ 2*(i3-1)+1 +  if(i4==3) 1 else 0 , "degree"] - sv_cn;
						SV_sum=SV_sum+sv_cn;
					}

				}else{


					 sv_cn=max(min (max(0, l_AS_CN[i1, "allele_1"] - l_AS_CN[i1 + if(i2==3) 1 else -1  , "allele_1"]), max(0, l_AS_CN[i3, "allele_1"] - l_AS_CN[i3 + if(i4==3) 1 else -1  , "allele_1"])) , min (max(0, l_AS_CN[i1, "allele_2"] - l_AS_CN[i1 + if(i2==3) 1 else -1  , "allele_2"]), max(0, l_AS_CN[i3, "allele_2"] - l_AS_CN[i3 + if(i4==3) 1 else -1  , "allele_2"])));
					 fi = which.max(c(min (max(0, l_AS_CN[i1, "allele_1"] - l_AS_CN[i1 + if(i2==3) 1 else -1  , "allele_1"]), max(0, l_AS_CN[i3, "allele_1"] - l_AS_CN[i3 + if(i4==3) 1 else -1  , "allele_1"])) , min (max(0, l_AS_CN[i1, "allele_2"] - l_AS_CN[i1 + if(i2==3) 1 else -1  , "allele_2"]), max(0, l_AS_CN[i3, "allele_2"] - l_AS_CN[i3 + if(i4==3) 1 else -1  , "allele_2"]))));
                                        if(!is.na(sv_cn)){
	                                       d[[fi]][ 2*(i1-1)+1 +  if(i2==3) 1 else 0 , "degree"] = d[[fi]][ 2*(i1-1)+1 +  if(i2==3) 1 else 0 , "degree"] - sv_cn;
	                                       d[[fi]][ 2*(i3-1)+1 +  if(i4==3) 1 else 0 , "degree"] = d[[fi]][ 2*(i3-1)+1 +  if(i4==3) 1 else 0 , "degree"] - sv_cn;
						SV_sum=SV_sum+sv_cn;
					}
				}
	#			print(l_AS_CN_degree[c(3,5),c(18,28,29,30)]);

			}
		}
		if(nrow(l_AS_CN) > 1){
			sum=0;
			for( d_i in 1:2){
				for (i in 1:((nrow(d[[d_i]]))/2 -1) ){
					j=i*2;
					if(!is.na(abs( d[[d_i]][j,3] - d[[d_i]][j+1,3]))){
						sum=sum+ abs( d[[d_i]][j,3] - d[[d_i]][j+1,3]);
					}
				}

			}
			acn_list = c(l_AS_CN$allele_1, l_AS_CN$allele_2);
			acn_list = unique(acn_list[!is.na(acn_list)]);
			min_acn_g_length=0;
			for(acn_list_i in  1:length(acn_list)){
				min_acn = acn_list[acn_list_i];
				min_acn_g_length = max(min_acn_g_length, sum(l_AS_CN[(l_AS_CN$allele_1 == min_acn | l_AS_CN$allele_2 == min_acn ), 3] - l_AS_CN[(l_AS_CN$allele_1 == min_acn | l_AS_CN$allele_2 == min_acn ), 2], na.rm=T));
			}
			print("here");
			print(sum(l_AS_CN$allele_1));
			print(sum(l_AS_CN$allele_2));
			print(l_AS_CN);
			print("sum");
			print(sum);
			print("SV_sum");
			print(SV_sum);
			print("g_length");
			print(min_acn_g_length)
	#		if(sum < diff || (sum == diff && SV_sum > SV_sum_max) || (sum == diff && SV_sum == SV_sum_max && min_acn_g_length  > g_length )) {

			edit_d_from_low_conf=0;
			edit_d_from_low_conf_switch=0;
			if(nrow(CN_low_conf)!=0){
	                        for( i in 1:nrow(CN_low_conf)){
					l_AS_CN_i=which(l_AS_CN[,2] == CN_low_conf[i,2]);
					edit1=abs(l_AS_CN[l_AS_CN_i,"allele_1"] -CN_low_conf[i,"allele_1"])+abs(l_AS_CN[l_AS_CN_i,"allele_2"] -CN_low_conf[i,"allele_2"]);
                                        edit2=abs(l_AS_CN[l_AS_CN_i,"allele_2"] -CN_low_conf[i,"allele_1"])+abs(l_AS_CN[l_AS_CN_i,"allele_1"] -CN_low_conf[i,"allele_2"]);
					if(!is.na(edit1) & !is.na(edit2)){
						edit_d_from_low_conf = edit_d_from_low_conf + min(edit1, edit2);
						if(edit1 != 0 && edit2 !=0 ){
							edit_d_from_low_conf_switch= edit_d_from_low_conf_switch+1;
						}
					}
				}
			}

			sum=sum+edit_d_from_low_conf_switch;
			print("switch_panalty_sum");
			print(sum);
	#		write.table(l_AS_CN,"l_AS_CN",quote=F, sep="\t", row.names=F);
	#		write.table(CN_low_conf, "CN_low_conf", quote=F, sep="\t", row.names=F);
		
#                        if(sum < r_diff || (sum == r_diff && SV_sum > SV_sum_max)) {
                        if(sum < r_diff ) {

				solution_count<<-1;
				r_diff <<- sum;
				SV_sum_max <<- SV_sum;
				AS_CN_opt <<- l_AS_CN;
				min_edit_d_from_low_conf <<- edit_d_from_low_conf;
				g_length <<- min_acn_g_length;
#			}else if (sum == r_diff && SV_sum == SV_sum_max){
                        }else if (sum == r_diff){

				solution_count<<-solution_count+1;
				if(edit_d_from_low_conf < min_edit_d_from_low_conf){
					AS_CN_opt <<- l_AS_CN;
					min_edit_d_from_low_conf <<- edit_d_from_low_conf;
					g_length <<- min_acn_g_length;
				}
			}
		}else{
			l_AS_CN$ACN_block=NA;
			if(!is.na(l_AS_CN[1,"allele_1"]) && l_AS_CN[1,"allele_1"] != l_AS_CN[1,"allele_2"]){
				l_AS_CN$ACN_block[1]=1;
			}
			AS_CN_opt <<- l_AS_CN;
		}
	       if(nrow(CN_low_conf) !=0){
			for( i in 1:nrow(CN_low_conf)){
				if(length(which(l_AS_CN[,2] == CN_low_conf[i,2] & l_AS_CN[,"allele_1"] == CN_low_conf[i,"allele_1"] & l_AS_CN[,"allele_2"] == CN_low_conf[i,"allele_2"]))!=0 ||
				  length(which(l_AS_CN[,2] == CN_low_conf[i,2] & l_AS_CN[,"allele_2"] == CN_low_conf[i,"allele_1"] & l_AS_CN[,"allele_1"] == CN_low_conf[i,"allele_2"]))!=0
				){
					CN_low_conf$solution_including_this[i] = CN_low_conf$solution_including_this[i] +1;
				}
			}
			CN_low_conf <<- CN_low_conf;
		}

		
	}

}

#save.image("1.RData");

save.image(file=paste(args[1],".RData",sep=""))

imp_multi_solution = imp[0,];
while(nrow(imp)>0){
	l_imp=imp[1,];
	imp=imp[-1,];
	prev_l_imp= l_imp[0,];
	while(nrow(l_imp) != nrow(prev_l_imp)  ){
		prev_l_imp= l_imp;
		if(nrow(imp)>0){
			imp_r=c();
			for(imp_i in 1:nrow(imp)){
				c1=c();
				c2=c();
				for(l_imp_i in 1:nrow(l_imp)){
					c1=c(c1,AS_CN[ max(1, l_imp[l_imp_i,1] -1 ):l_imp[l_imp_i,2] , 3]);
					c2=c(c2,AS_CN[ l_imp[l_imp_i,1] : min ( nrow(AS_CN), l_imp[l_imp_i,2] +1 ) , 2]);
				}
				c3=AS_CN[ max(1, imp[imp_i,1] -1 ) : imp[imp_i,2], 3];
				c4=AS_CN[ imp[imp_i,1] : min( nrow(AS_CN), imp[imp_i,2] + 1 ), 2];
	
				exist=length(which((SV[,3] %in% c(c1,c2) & SV[,5] %in% c(c3,c4)) |(SV[,3] %in% c(c3,c4) & SV[,5] %in% c(c1,c2) )));
				if(exist){
					l_imp = rbind (l_imp, imp[imp_i,]);
					imp_r=c(imp_r, imp_i);
				}
			}
			if(length(imp_r)>0){
				imp=imp[-imp_r,];
			}
		}
	}
	solution_count=0;
	r_diff=Inf;
	SV_sum_max=0;
	g_length=0;
	AS_CN_opt=AS_CN;
	search_space = 1;
       for(search_i in 1:nrow(l_imp)){
		search_space=search_space*length(which(!is.na(l_imp[search_i, 3:ncol(l_imp)])))
	}
	p_thres_h = p_thres;
	while(search_space > imp_cut_thres && p_thres_h > 0){
		p_thres_h = p_thres_h -0.05;
		for(search_i in 1:nrow(l_imp)){
			temp=CN_low_conf[CN_low_conf$Start.bp %in% AS_CN[l_imp[search_i,1]:l_imp[search_i,2],]$Start.bp,]
			max_allele_i=3;
			max_allele=0;
			for(search_j in 3:ncol(l_imp)){
				if(max_allele <= length(unique(which(l_imp[search_i,search_j] == temp[temp$p > p_thres_h,"allele_1"]), which(l_imp[search_i,search_j] == temp[temp$p > p_thres_h,"allele_2"])))){
				max_allele_i=search_j;
				max_allele= length(unique(which(l_imp[search_i,search_j] == temp[temp$p > p_thres_h,"allele_1"]), which(l_imp[search_i,search_j] == temp[temp$p > p_thres_h,"allele_2"])));
				}
			}
			if(max_allele !=0){
	                        for(search_j in 3:ncol(l_imp)){
					if(search_j !=max_allele_i){
						l_imp[search_i, search_j] = NA;
					}
				}
			}
		}
		
		search_space=1;
		for(search_i in 1:nrow(l_imp)){
			search_space=search_space*length(which(!is.na(l_imp[search_i, 3:ncol(l_imp)])))
		}
	}

	recursive_imp (AS_CN, 1, nrow(l_imp), l_imp);
	AS_CN=AS_CN_opt;


	if(solution_count > 1){
		imp_multi_solution = rbind(imp_multi_solution, l_imp);
	}
}


r_diff=Inf;
SV_sum_max=0;
g_length=0;
AS_CN_opt=AS_CN;
solution_count=0;
min_edit_d_from_low_conf=Inf;
recursive_imp (AS_CN, Inf, -Inf, imp);
AS_CN=AS_CN_opt;




rm(i)
for( CN_i in 1:nrow(AS_CN)){
        if(1){
#	if(AS_CN[CN_i,"allele_1"] !=0 && AS_CN[CN_i,"allele_2"] !=0 && is.na(AS_CN[CN_i,"mean"])){
                t=SNP[SNP$V1==AS_CN[CN_i,1]&SNP$V2>AS_CN[CN_i,2]&SNP$V2<AS_CN[CN_i,3],];
                t=t[t$V8< MAX_BASE_DEPTH & t$V9< MAX_BASE_DEPTH,];
                if(nrow(t)!=0){
                        t_IQR=IQR(c(t[,8],t[,9]));
                        Q1=quantile(c(t[,8],t[,9]), 1/4);
                        Q3=quantile(c(t[,8],t[,9]), 3/4);
                        t=t[t$V8 > Q1-coeff_IQR*t_IQR & t$V8 < Q3+coeff_IQR*t_IQR & t$V9 > Q1-coeff_IQR*t_IQR & t$V9 < Q3+coeff_IQR*t_IQR,]
                }

		if(nrow(t)> 1){
                        estimate_hat=fitdist(t$V8+t$V9, "nbinom")$estimate
                        NB_MEAN=estimate_hat[2]/AS_CN[CN_i,18];
                        NB_DIS=estimate_hat[1];
                        ML_ACN=NA;
                        ML_sum=-Inf;
			ACN=min(AS_CN[CN_i,"allele_1"] , AS_CN[CN_i,"allele_2"]);

			if(AS_CN[CN_i,"allele_1"] != AS_CN[CN_i,"allele_2"]){

                                        NB_mean=NB_MEAN;
                                        NB_dis1=NB_DIS;
					NB_dis2=NB_DIS;
                                        NB_value=Inf;
                                        E=matrix(ncol=4, nrow=nrow(t))
                                        max_iter=MAX_ITER;
                                        NB_value_diff=Inf;
                                        while(max_iter ==MAX_ITER || ( max_iter > 0 && NB_value_diff > EM_tol*NB_value)) {
                                                total_p=dnbinom(t[,8], size=NB_dis1, mu=(ACN*purity+1*(1-purity))*NB_mean)*dnbinom(t[,9], size=NB_dis2, mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*NB_mean)+
                                                        dnbinom(t[,8], size=NB_dis2, mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*NB_mean)*dnbinom(t[,9], size=NB_dis1, mu=(ACN*purity+1*(1-purity))*NB_mean);

                                                E[,1] = dnbinom(t[,8], size=NB_dis1, mu=(ACN*purity+1*(1-purity))*NB_mean)*dnbinom(t[,9], size=NB_dis2, mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*NB_mean)/total_p;
                                                E[,2] = 1-E[,1];


                                                frr<-function(x){
                                                        sum=0;
                                                        sum=sum(E[,1]*log(dnbinom(t[,8],mu=(ACN*purity+1*(1-purity))*x[1], size=x[2])*dnbinom(t[,9], mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*x[1], size=x[3]))+
                                                                E[,2]*log(dnbinom(t[,8], mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*x[1], size=x[3])*dnbinom(t[,9],mu=(ACN*purity+1*(1-purity))*x[1], size=x[2])), na.rm=T)
                                                        -sum;
                                                }

                                                optim_err=0;
                                                rm(M);


                                        tryCatch({M=optim(c(NB_MEAN,NB_DIS,NB_DIS),frr, method="L-BFGS-B", control=list(factr=optim_tol), lower=1)
                                                }, error=function(e){optim_err=1;});

                                                if(optim_err==0 && exists("M")){
                                                        print(M);
                                                        NB_value_diff=abs(NB_value-M$value);
                                                        NB_mean=M$par[1];
                                                        NB_dis1=M$par[2];
							NB_dis2=M$par[3];
                                                        NB_value=M$value;
##                                                      NB_value = NB_frr(M$par);
                                                        max_iter=max_iter-1;
                                                }else{
                                                        max_iter=-1;
                                                }

                                        }

                                        mixture_L_sum=0;
                                        for(m_i in 1:nrow(t)){
                                                mixture_L_sum = mixture_L_sum + max (dnbinom(t[m_i,8],mu=(ACN*purity+1*(1-purity))*NB_mean, size=NB_dis1)*dnbinom(t[m_i,9], mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*NB_mean, size=NB_dis2),dnbinom(t[m_i,9],mu=(ACN*purity+1*(1-purity))*NB_mean, size=NB_dis1)*dnbinom(t[m_i,8], mu=((CN[CN_i,18]-ACN)*purity+1*(1-purity))*NB_mean, size=NB_dis2));
                                        }


                                        if( !is.na(mixture_L_sum)){
						AS_CN[CN_i, "mean"]=NB_mean;
						AS_CN[CN_i,"dis1"]=NB_dis1;
						AS_CN[CN_i,"dis2"]=NB_dis2;
                                         }


			}else{
                                        NB_mean=NB_MEAN;
                                        NB_dis=NB_DIS;

                                        frr<-function(x){
                                                sum=0;
                                                sum=sum(log(dnbinom(t[,8],mu=(ACN*purity+1*(1-purity))*x[1], size=x[2])*dnbinom(t[,9],mu=(ACN*purity+1*(1-purity))*x[1], size=x[2])), na.rm=T);
                                                -sum;
                                       }

                                        optim_err=0;
                                        rm(M);
                                        tryCatch({M=optim(c(NB_MEAN,NB_DIS),frr, method="L-BFGS-B", control=list(factr=optim_tol), lower=1)
                                                }, error=function(e){optim_err=1;});

                                        if(optim_err==0 && exists("M")){
                                                print(M);
                                                NB_value_diff=abs(NB_value-M$value);
                                                NB_mean=M$par[1]
                                                NB_dis=M$par[2]
                                                NB_value=M$value;
                                                max_iter=max_iter-1;
                                        }else{
                                                max_iter=-1;
                                        }
                                        mixture_L_sum = 0;
                                        mixture_L_sum=sum((dnbinom(t[,8],mu=(ACN*purity+1*(1-purity))*NB_mean, size=NB_dis)*dnbinom(t[,9],mu=(ACN*purity+1*(1-purity))*NB_mean, size=NB_dis)));

                                        if( !is.na(mixture_L_sum)){
                                                AS_CN[CN_i, "mean"]=NB_mean;
                                                AS_CN[CN_i,"dis1"]=NB_dis;
						AS_CN[CN_i,"dis2"]=NB_dis;
                                         }


			}
		}
	}
}


write.table(AS_CN, paste("copy_numbers.CN_opt.ACN.",args[1],sep=""), sep="\t", quote=F, row.names=F)
