args <- commandArgs(TRUE)
library(lpSolve)
library(lpSolveAPI)
####
#args=c();
#args[1]="simul"
#marker_size=2;
marker_size=0;
false_marker_size=1;
min_segment_length=50;
ABSOLUTE_subclonal_ccf_ci95_high=0.8;
min_expected_modal_cn_diff=0.3;
#ABSOLUTE_CN_limit=14;
#ABSOLUTE_copy_ratio_limit=2.5;
static_CN_limit=Inf;
read_length=args[4];
marker_sparsity=0;
#0#ABSOLUTE_CN_limit=Inf; ####################### nearest regression is performed
modal_cn_delta=1;
modal_cn_delta_for_high_copy=0;
modal_cn_delta_for_zero_marker=0;
lp_obj_thres=2;
lp_iter_thres=100;
##################  x <  impute_lv + impute_lv_plus
impute_lv=2; 
impute_lv_plus=2;
max_SV_min=5000;
time_out=240;
SV_max=200; ######## upper bound of a single SV
CN_max=200; ####### for initial relaxed  problem
seg.q.tab.threshold=0.95
####
load(paste("ABSOLUTE_output/",args[1],".ABSOLUTE.RData",sep=""));
ABS_purity_ploidy=read.table(paste("ABSOLUTE_output/output/reviewed/",args[1],".test.ABSOLUTE.table.txt",sep=""), header=T,sep="\t");
ABS_purity=ABS_purity_ploidy$purity;
ABS_ploidy=ABS_purity_ploidy$ploidy;
ABS_D=ABS_purity*ABS_ploidy+2*(1-ABS_purity);

ABS_output=read.table(paste("ABSOLUTE_output/output/reviewed/SEG_MAF/",args[1],".segtab.txt",sep=""), header=T)
global_copy_ratio=sum(ABS_output$W*ABS_output$copy_ratio,na.rm=T)*2;
global_Log=log(global_copy_ratio)/log(2);

ABS_output$seg.q.tab=0;
for(i in 1:nrow(ABS_output)){
	ABS_output$seg.q.tab[i]=max(seg.dat$mode.res$seg.q.tab[1,i,1:15]);
}
##ABS_output$seg.q.tab=seg.dat$mode.res$seg.q.tab;
ABS_output_duplicated=ABS_output;
#ABS_output[ABS_output$expected_cn>=ABSOLUTE_CN_limit,"expected_cn"]=Inf;
ABS_input=read.table(args[2], header=T)
colnames(ABS_input)=c("Chromosome","Start.bp","End.bp", "Probes", "Segment_Mean")
new=merge(ABS_input, ABS_output, by=c("Chromosome","Start.bp", "End.bp"), all.x=T)

new=new[order(new[,1],new[,2]),];

new$raw_expected_cn=new$expected_cn;
new$Segment_Mean=new$Segment_Mean-global_Log;
ABS_input$Segment_Mean=ABS_input$Segment_Mean-global_Log;


for(i in 1:nrow(new)){
	copy_ratio=(2^ABS_input[ABS_input$Chromosome==new$Chromosome[i] & ABS_input$Start==new$Start.bp[i] & ABS_input$End==new$End.bp[i],"Segment_Mean"]);
        int_CN=max(0,2*(copy_ratio*ABS_D/2/ABS_purity-(1-ABS_purity)/ABS_purity));

	new[i,"raw_expected_cn"]=int_CN;
	new[i,"expected_cn"]=int_CN;
}

new$alternative_modal_cn=new$modal_cn;
#new$Probes=new$n_probes;


for(i in 1:nrow(new)){
#        if(is.na(new$modal_cn[i]) || new$End.bp[i]-new$Start.bp[i] < min_segment_length || new$Probes[i]<=0){
        if(is.na(new$Segment_Mean[i]) || new$End.bp[i]-new$Start.bp[i] < min_segment_length || new$Probes[i]<=0){

                new$modal_cn[i]=NA;
                new[i,"expected_cn"]=Inf;
#        }else if(new$copy_ratio[i]==ABSOLUTE_copy_ratio_limit){
        }else if(abs(new$modal_cn[i]- new$raw_expected_cn[i]) > 1 ){

                copy_ratio=(2^ABS_input[ABS_input$Chromosome==new$Chromosome[i] & ABS_input$Start==new$Start.bp[i] & ABS_input$End==new$End.bp[i],"Segment_Mean"]);
                int_CN=max(0,2*(copy_ratio*ABS_D/2/ABS_purity-(1-ABS_purity)/ABS_purity));

                new[i,"raw_expected_cn"]=int_CN;##############################
         #       new[i,"expected_cn"]=Inf;
###############################################################################
                new[i,"modal_cn"]=min(CN_max,round(int_CN));
	}else if(new[i,"n_probes"] < 0){
    #            new[i,"expected_cn"]=Inf;
        }else if(new[i,"n_probes"]<=marker_size){
                new[i,"expected_cn"]=Inf;
        }else if(new[i,"seg.q.tab"]<seg.q.tab.threshold){
                new[i,"expected_cn"]=Inf;
#                copy_ratio=2*(2^ABS_input[ABS_input$Chromosome==new$Chromosome[i] & ABS_input$Start==new$Start.bp[i] & ABS_input$End==new$End.bp[i],"Segment_Mean"]);
#                int_CN=copy_ratio*ABS_D/2/ABS_purity-(1-ABS_purity)/ABS_purity;
#                new[i,"raw_expected_cn"]=int_CN;
        }else{
#                copy_ratio=2*(2^ABS_input[ABS_input$Chromosome==new$Chromosome[i] & ABS_input$Start==new$Start.bp[i] & ABS_input$End==new$End.bp[i],"Segment_Mean"]);
#                int_CN=copy_ratio*ABS_D/2/ABS_purity-(1-ABS_purity)/ABS_purity;
#                new[i,"expected_cn"]=int_CN;
        }
}

for(i in 1:nrow(new)){
        if(!is.na(new[i,"raw_expected_cn"]) && !is.na(new[i,"modal_cn"])){
                alternatives=c(floor(new[i,"raw_expected_cn"]),ceiling(new[i,"raw_expected_cn"]));
		alternative_modal_cns=alternatives[which(alternatives!=new[i,"modal_cn"])];
		if(length(alternative_modal_cns) > 0 ){
	                new$alternative_modal_cn[i]=alternatives[which(alternatives!=new[i,"modal_cn"])][1];
		}else{
                        new$alternative_modal_cn[i]=new$modal_cn[i];
		}
                if(abs(new[i,"raw_expected_cn"]-new[i,"modal_cn"])>min_expected_modal_cn_diff){
                       new[i,"seg.q.tab"]=0;
	                new[i,"expected_cn"]=Inf;
                }
                if(abs(new[i,"Probes"])*as.integer(read_length)/abs(new[i,"End.bp"]-new[i,"Start.bp"])<marker_sparsity){
                        new[i,"seg.q.tab"]=0;
	                new[i,"expected_cn"]=Inf;
                }
        }
}


write.table(new,"new_debug",quote=F, sep="\t", row.names=F)
SV_table=read.table(args[3])
levels(SV_table[,2])[levels(SV_table[,2])=="X"]<-"23"
levels(SV_table[,4])[levels(SV_table[,4])=="X"]<-"23"
SV_table[,2]=as.integer(as.character(SV_table[,2]))
SV_table[,4]=as.integer(as.character(SV_table[,4]))
#SV_table$V2 <- factor(SV_table$V2, levels = paste(sort(as.integer(levels(SV_table$V2)))))
#SV_table$V4 <- factor(SV_table$V4, levels = paste(sort(as.integer(levels(SV_table$V4)))))



SV_set=SV_table;

 SV_matrix_minus <- function(A,B){
	if(nrow(B)==0)
		return(A)
	test=c()
	for(i in 1:nrow(A)){
		for(j in 1:nrow(B)){
			if(A[i,2]==B[j,2]&&A[i,3]==B[j,3]&&A[i,4]==B[j,4]&&A[i,5]==B[j,5])
				test=c(test,i)
		}
	}
	return(A[-test,])
}

SV_constraint_explicit <-function (con_A, con_A_variable, con_solution){
        con_not_satisfied=0;
        for(con_A_i in 1:nrow(con_A)){
                if(con_A_variable[con_A_i,4]==0){
                        con_A_4=con_A[con_A_i,4];
                }else{
                        con_A_4=con_solution[con_A_variable[con_A_i,4]];
                }
                if(con_A_variable[con_A_i,5]==0){
                        con_A_5=con_A[con_A_i,5];
                }else{
                        con_A_5=con_solution[con_A_variable[con_A_i,5]];
                }
                if(con_A_variable[con_A_i,6]==0){
                        con_A_6=con_A[con_A_i,6];
                }else{
                        con_A_6=con_solution[con_A_variable[con_A_i,6]];
                }
                if(con_A_variable[con_A_i,7]==0){
                        con_A_7=con_A[con_A_i,7];
                }else{
                        con_A_7=con_solution[con_A_variable[con_A_i,7]];
                }
                if(con_A_5>max(con_A_4-con_A_6,0) || con_A_7>max(con_A_6-con_A_4,0))
                        con_not_satisfied=con_not_satisfied+1;

        }
        if(con_not_satisfied!=0){
                return(FALSE)
        }else{
                return(TRUE);
        }
 }


 not_in_queue <- function(e1, queue){
	if(length(queue)==0)
		return(T);
	for(queue_i in 1:length(queue)){
		if(all(e1==queue[[queue_i]])==T)
			return(F)
	}
	return(T);
 }
 add_A <- function(A,node){
        right_search=new[new$Chromosome==node[1]&new$Start.bp==node[2]+1,];
        left_search=new[new$Chromosome==node[1]&new$End.bp==node[2]-1,];
        if(nrow(right_search)!=0){
                A=rbind(A,c(right_search[1,1],node[2],right_search[1,2],0,0,0,0))
        }else if(nrow(left_search)!=0){
                A=rbind(A,c(left_search[1,1],left_search[1,3],node[2],0,0,0,0))
        }
	return(A);


 }


 related_set_find_BFS_ver <- function (SV_line, SV, A){
	
	A_SV_search=c();
	SV_duplicate=SV;
	BFS_queue=list();
	Empty_set=list()
        s=strsplit(as.character(SV_line[1,6]),"to");
	BFS_root=c(as.integer(as.character(SV_line[1,2])),SV_line[1,3]);
	BFS_queue[[1]]=c(as.integer(as.character(SV_line[1,2])),SV_line[1,3]);
	Empty_set[[1]]=c(as.integer(as.character(SV_line[1,2])),SV_line[1,3]);
	A=add_A(A,BFS_root);

	while(length(BFS_queue)!=0){
		BFS_current=BFS_queue[[1]];
		BFS_queue=BFS_queue[-1];
		################ search for neighbor
		BFS_current_neighbor=list();
		Start_search=new[new$Chromosome==BFS_current[1]&new$Start.bp==BFS_current[2],];
		if(nrow(Start_search)!=0){
			current_left=c(new[new$Chromosome==BFS_current[1] & new$End.bp==BFS_current[2]-1,1], new[new$Chromosome==BFS_current[1] & new$End.bp==BFS_current[2]-1,3]);
			if(length(current_left)!=0)
				BFS_current_neighbor[[length(BFS_current_neighbor)+1]]=current_left;
			if(is.na(Start_search$Segment_Mean[1])|| Start_search$Probes[1]<=marker_size || Start_search$expected_cn[1] == Inf || abs(Start_search$expected_cn[1]-Start_search$modal_cn[1]) > min_expected_modal_cn_diff)
				BFS_current_neighbor[[length(BFS_current_neighbor)+1]]=c(Start_search[1,1], Start_search[1,3])
		}
                End_search=new[new$Chromosome==BFS_current[1]&new$End.bp==BFS_current[2],];
                if(nrow(End_search)!=0){
                        current_right=c(new[new$Chromosome==BFS_current[1] & new$Start.bp==BFS_current[2]+1,1], new[new$Chromosome==BFS_current[1] & new$Start.bp==BFS_current[2]+1,2]);
			if(length(current_right)!=0)
                                BFS_current_neighbor[[length(BFS_current_neighbor)+1]]=current_right;
                        if(is.na(End_search$Segment_Mean[1])|| End_search$Probes[1]<=marker_size || End_search$expected_cn[1] == Inf ||abs(End_search$expected_cn[1]-End_search$modal_cn[1]) > min_expected_modal_cn_diff)
                                BFS_current_neighbor[[length(BFS_current_neighbor)+1]]=c(End_search[1,1], End_search[1,2])
		}
		which_SV_V3_search=which(SV$V2==BFS_current[1]&SV$V3==BFS_current[2]);
                duplicate_which_SV_V3_search=which(SV_duplicate$V2==BFS_current[1]&SV_duplicate$V3==BFS_current[2]);
		SV_V3_search=SV[which_SV_V3_search,];
		if(nrow(SV_V3_search)!=0){
			SV=SV[-which_SV_V3_search,];
			A_SV_search=c(A_SV_search,duplicate_which_SV_V3_search);
			for(SV_V3_search_i in 1:nrow(SV_V3_search)){
				BFS_current_neighbor[[length(BFS_current_neighbor)+1]]=c(SV_V3_search[SV_V3_search_i,4],SV_V3_search[SV_V3_search_i,5]);
			}
		}
		which_SV_V5_search=which(SV$V4==BFS_current[1]&SV$V5==BFS_current[2]);
                duplicate_which_SV_V5_search=which(SV_duplicate$V4==BFS_current[1]&SV_duplicate$V5==BFS_current[2]);
		SV_V5_search=SV[which_SV_V5_search,];
		if(nrow(SV_V5_search)!=0){
                        SV=SV[-which_SV_V5_search,];
                        A_SV_search=c(A_SV_search,duplicate_which_SV_V5_search);
                        for(SV_V5_search_i in 1:nrow(SV_V5_search)){
                        	BFS_current_neighbor[[length(BFS_current_neighbor)+1]]=c(SV_V5_search[SV_V5_search_i,2],SV_V5_search[SV_V5_search_i,3]);
			}
		}	
		#################
		if(length(BFS_current_neighbor)!=0){
			for(neighbor_i in 1:length(BFS_current_neighbor)){
				if(not_in_queue(BFS_current_neighbor[[neighbor_i]],Empty_set)){
				        Empty_set[[length(Empty_set)+1]]=BFS_current_neighbor[[neighbor_i]];
					A=add_A(A,BFS_current_neighbor[[neighbor_i]]);
					BFS_queue[[length(BFS_queue)+1]]=BFS_current_neighbor[[neighbor_i]];
				}
					
			}
		}

	}
        A=A[order(A[,1],A[,2]),]
	A=unique(A);
        return(list(A,SV,A_SV_search));


 } 


 related_set_find <- function (SV_line,SV,A){
 	s=strsplit(as.character(SV_line[1,6]),"to");
	if(s[[1]][1]==3){
		A=rbind(A,c(as.integer(as.character(SV_line[1,2])),SV_line[1,3],SV_line[1,3]+1, 0, 0, 0, 0));
		}else
		A=rbind(A, c(as.integer(as.character(SV_line[1,2])), SV_line[1,3]-1, SV_line[1,3], 0, 0, 0, 0));
	if(s[[1]][2]==3){
                A=rbind(A,c(as.integer(as.character(SV_line[1,4])),SV_line[1,5],SV_line[1,5]+1, 0, 0, 0, 0));
	}else
                A=rbind(A, c(as.integer(as.character(SV_line[1,4])), SV_line[1,5]-1, SV_line[1,5], 0, 0, 0, 0));

	prev_nrow_A=0;
	while(prev_nrow_A!=nrow(A)){
		prev_nrow_A=nrow(A);
		for(j in 1:nrow(A)){
			left_index=TRUE;
			left_pos=A[j,2];
			while(left_index==T){

## marker size < 100 ########
				##print(left_pos);
				if(nrow(new[new$End.bp==left_pos,])!=0&&nrow(new[new$End.bp==(new[new$End.bp==left_pos,][1,2]-1),])!=0&&(is.na(new[new$End.bp==left_pos,5])==T||new[new$End.bp==left_pos,4]<=marker_size||(new[new$End.bp==left_pos,20]==1&&new[new$End.bp==left_pos,23]>ABSOLUTE_subclonal_ccf_ci95_high))){
					A=rbind(A,c(as.integer(as.character(new[new$End.bp==left_pos,][1,1])),new[new$End.bp==left_pos,][1,2]-1, new[new$End.bp==left_pos,][1,2], 0, 0, 0, 0));
					left_pos=new[new$End.bp==left_pos,][1,2]-1;
				}else{
					left_index=FALSE;
				}
			}
               		right_index=TRUE;
                	right_pos=A[j,3];
                	while(right_index==T){
                        	if(nrow(new[new$Start.bp==right_pos,])!=0&&nrow(new[new$Start.bp==(new[new$Start.bp==right_pos,][1,3]+1),])!=0&&(is.na(new[new$Start.bp==right_pos,5])==T||new[new$Start.bp==right_pos,4]<=marker_size||(new[new$Start.bp==right_pos,20]==1&&new[new$Start.bp==right_pos,23]>ABSOLUTE_subclonal_ccf_ci95_high))){
                                	A=rbind(A,c(as.integer(as.character(new[new$End.bp==left_pos,][1,1])),new[new$Start.bp==right_pos,][1,3], new[new$Start.bp==right_pos,][1,3]+1, 0, 0, 0, 0));
                                	right_pos=new[new$Start.bp==right_pos,][1,3]+1;
                        	}else{
                                	right_index=FALSE;
                        	}
                	}
		}
		A=unique(A);
	}

	A=unique(A);

        prev_nrow_A=0;
        while(prev_nrow_A!=nrow(A)){
	        prev_nrow_A=nrow(A);

		for(j in 1:nrow(A)){
			t=SV[SV$V2==A[j,1]&SV$V3==A[j,2],];
			SV<-SV[SV$V2!=A[j,1]|SV$V3!=A[j,2],];
			if(is.na(t[1,1])!=T){
				for(k in 1:nrow(t)){
                                        related_set_find_return=related_set_find(t[k,],SV,A);
                                        A=related_set_find_return[[1]];
                                        SV=related_set_find_return[[2]];
				}
			}

                	t=SV[SV$V4==A[j,1]&SV$V5==A[j,2],];
			SV<-SV[SV$V4!=A[j,1]|SV$V5!=A[j,2],];
               		if(is.na(t[1,1])!=T){
               	        	for(k in 1:nrow(t)){
                                        related_set_find_return=related_set_find(t[k,],SV,A);
                                        A=related_set_find_return[[1]];
                                        SV=related_set_find_return[[2]];

				 }
               		}
                	t=SV[SV$V2==A[j,1]&SV$V3==A[j,3],];
                	SV<-SV[SV$V2!=A[j,1]|SV$V3!=A[j,3],];
                	if(is.na(t[1,1])!=T){
                        	for(k in 1:nrow(t)){
                                        related_set_find_return=related_set_find(t[k,],SV,A);
                                        A=related_set_find_return[[1]];
                                        SV=related_set_find_return[[2]];
                        	}
                	}
                	t=SV[SV$V4==A[j,1]&SV$V5==A[j,3],];
                	SV<-SV[SV$V4!=A[j,1]|SV$V5!=A[j,3],];
                	if(is.na(t[1,1])!=T){
                        	for(k in 1:nrow(t)){
                                        related_set_find_return=related_set_find(t[k,],SV,A);
                                        A=related_set_find_return[[1]];
                                        SV=related_set_find_return[[2]];
                       		}
                	}

		}
	}
	return(list(unique(A),SV))

}

Integer_prog_H_positive <- function(L, nrow_A_CN){
        local_l=L[[1]];
        local_A=L[[2]];
	f.obj=c();
        A_variable=local_A;
        A_variable[1:nrow(local_A),4:7]=0;
        for(l_i in 1:length(local_l)){
                f.obj[l_i]=0;
                for(l_i_2 in 1:ncol(local_l[[l_i]][[2]])){
                        if(local_l[[l_i]][[2]][2,l_i_2]== 4 || local_l[[l_i]][[2]][2,l_i_2]== 6){
                                f.obj[l_i]=f.obj[l_i]+1;
                        }else{
                                f.obj[l_i]=f.obj[l_i]-1;
                        }
                        A_variable[local_l[[l_i]][[2]][1,l_i_2],local_l[[l_i]][[2]][2,l_i_2]]=l_i;
                }
        }
        for(l_i in (length(local_l)+1):(length(local_l)+nrow(local_A))){
                f.obj[l_i]=-2;
        }
	modal_cn_both_bound_count=0;
	modal_cn_lower_bound_count_for_absolute_max_cn=0;
	if(nrow_A_CN!=0){
	        for(l_i in 1:nrow_A_CN){
##	                if(!is.na(l[[l_i]][[5]][4]) && l[[l_i]][[5]][4]>=ABSOLUTE_CN_limit){
                        if(!is.na(l[[l_i]][[5]][4]) && l[[l_i]][[5]][4]>=Inf){
				modal_cn_lower_bound_count_for_absolute_max_cn=modal_cn_lower_bound_count_for_absolute_max_cn+1;
			}else if(!is.na(l[[l_i]][[5]][4])){
				modal_cn_both_bound_count=modal_cn_both_bound_count+2;
			}
	        }
	}

#	f.con=matrix(0, nrow=nrow(A)+length(l), ncol=length(l));
#####################
	precount=0;
	SV_const_count=0;
	SV_SOS_count=0;
#	if(nrow_A_CN!=0){
#		for(l_i in 1:nrow_A_CN){
#                	for(l_i_c in 1:ncol(local_l[[l_i]][[2]])){
#				precount=precount+4;
#			}
#		}
#	}else{
#		precount=0;
#	}
#############################
for(A_variable_i in 1:nrow(A_variable)){
	if(A_variable[A_variable_i,5]!=0){
		if(!(A_variable[A_variable_i,4]==0 && A_variable[A_variable_i,6]==0 && A_variable[A_variable_i,7]==0)){
			SV_const_count=SV_const_count+2;
			SV_SOS_count=SV_SOS_count+2;
		}else{
			SV_const_count=SV_const_count+1;
		}
	}
        if(A_variable[A_variable_i,7]!=0){
                if(!(A_variable[A_variable_i,4]==0 && A_variable[A_variable_i,6]==0 && A_variable[A_variable_i,5]==0)){
                        SV_const_count=SV_const_count+2;
                        SV_SOS_count=SV_SOS_count+2;
                }else{
                        SV_const_count=SV_const_count+1;
                }
	}
}
###############################
        f.con=matrix(0, nrow=nrow(local_A)*2+nrow_A_CN+precount+(length(local_l)-nrow_A_CN)+SV_const_count+modal_cn_both_bound_count+modal_cn_lower_bound_count_for_absolute_max_cn+1, ncol=length(local_l)+nrow(local_A)+SV_SOS_count);
	f.con_i=1;
	f.dir=c();
	f.dir[1:(nrow(local_A)*2)]=rep(">=", (nrow(local_A)*2));
	if(nrow_A_CN!=0)
		f.dir[(nrow(local_A)*2+1):(nrow(local_A)*2+nrow_A_CN)]=rep("<=", nrow_A_CN);
#	if(nrow_A_CN!=0)
#		f.dir[(nrow(local_A)*2+nrow_A_CN+1):(nrow(local_A)*2+nrow_A_CN+precount)]=rep("<", precount);
	f.dir[(nrow(local_A)*2+nrow_A_CN+1):(nrow(local_A)*2+nrow_A_CN+(length(local_l)-nrow_A_CN))]=rep("<=", length(local_l)-nrow_A_CN);

	if(SV_const_count!=0)
		f.dir[(nrow(local_A)*2+nrow_A_CN+(length(local_l)-nrow_A_CN)+1):(nrow(local_A)*2+nrow_A_CN+(length(local_l)-nrow_A_CN)+SV_const_count)]=rep("<=", SV_const_count);

	if(modal_cn_both_bound_count!=0)
		f.dir[(nrow(local_A)*2+nrow_A_CN+precount+(length(local_l)-nrow_A_CN)+SV_const_count+1):(nrow(local_A)*2+nrow_A_CN+precount+(length(local_l)-nrow_A_CN)+SV_const_count+modal_cn_both_bound_count)]=rep(c(">=","<="),modal_cn_both_bound_count/2);

	if(modal_cn_lower_bound_count_for_absolute_max_cn!=0)
		f.dir[(nrow(local_A)*2+nrow_A_CN+precount+(length(local_l)-nrow_A_CN)+SV_const_count+modal_cn_both_bound_count+1):(nrow(local_A)*2+nrow_A_CN+precount+(length(local_l)-nrow_A_CN)+SV_const_count+modal_cn_both_bound_count+modal_cn_lower_bound_count_for_absolute_max_cn)]=rep(">=",modal_cn_lower_bound_count_for_absolute_max_cn);

	f.dir[nrow(f.con)]=">=";
	f.rhs=c();
	f.rhs_i=1;

	if(SV_SOS_count!=0)
		f.obj[(length(local_l)+nrow(local_A)+1):(length(local_l)+nrow(local_A)+SV_SOS_count)]=0;

	for(A_i in 1:nrow(A_variable)){
		f.rhs[f.rhs_i]=0;
		if(A_variable[A_i,4]==0){
			f.rhs[f.rhs_i]=f.rhs[f.rhs_i]-local_A[A_i,4];
		}else{
			f.con[f.con_i,A_variable[A_i,4]]=f.con[f.con_i,A_variable[A_i,4]]+1;
		}
                if(A_variable[A_i,5]==0){
                        f.rhs[f.rhs_i]=f.rhs[f.rhs_i]-local_A[A_i,5];
                }else{
                        f.con[f.con_i,A_variable[A_i,5]]=f.con[f.con_i,A_variable[A_i,5]]-1;
                }
                f.con[f.con_i,length(local_l)+A_i]=-1;
                f.rhs_i=f.rhs_i+1;
                f.con_i=f.con_i+1;

                f.rhs[f.rhs_i]=0;
		if(A_variable[A_i,6]==0){
                        f.rhs[f.rhs_i]=f.rhs[f.rhs_i]-local_A[A_i,6];
                }else{
                        f.con[f.con_i,A_variable[A_i,6]]=f.con[f.con_i,A_variable[A_i,6]]+1;
                }
                if(A_variable[A_i,7]==0){
                        f.rhs[f.rhs_i]=f.rhs[f.rhs_i]-local_A[A_i,7];
                }else{
                        f.con[f.con_i,A_variable[A_i,7]]=f.con[f.con_i,A_variable[A_i,7]]-1;
                }

		f.con[f.con_i,length(local_l)+A_i]=-1;
		f.rhs_i=f.rhs_i+1;
		f.con_i=f.con_i+1;
	}
	if(nrow_A_CN!=0){
		for(l_i in 1:nrow_A_CN){
	################ lower bound = 1 or 0
	#                f.rhs[f.rhs_i]=1;
			f.rhs[f.rhs_i]=CN_max;
	################
			f.con[f.con_i,l_i]=1;
			f.rhs_i=f.rhs_i+1;
			f.con_i=f.con_i+1;
		}

	}
############ upper bound of single SV #########
	for(l_i in (nrow_A_CN+1):length(local_l)){
		f.rhs[f.rhs_i]=SV_max;
		f.con[f.con_i,l_i]=1;
                f.rhs_i=f.rhs_i+1;
                f.con_i=f.con_i+1;

	}
###########################################a
SOS_i=1;
for(A_variable_i in 1:nrow(A_variable)){
        if(A_variable[A_variable_i,5]!=0){
                if(!(A_variable[A_variable_i,4]==0 && A_variable[A_variable_i,6]==0 && A_variable[A_variable_i,7]==0)){
			f.con[f.con_i,A_variable[A_variable_i,5]]=+1;
			f.rhs[f.rhs_i]=0;
			if(A_variable[A_variable_i,4]==0){
				f.rhs[f.rhs_i]=f.rhs[f.rhs_i]+local_A[A_variable_i,4];
			}else{
				f.con[f.con_i,A_variable[A_variable_i,4]]=-1;
			}
                        if(A_variable[A_variable_i,6]==0){
                                f.rhs[f.rhs_i]=f.rhs[f.rhs_i]-local_A[A_variable_i,6];
                        }else{
                                f.con[f.con_i,A_variable[A_variable_i,6]]=+1;
                        }
			if(A_variable[A_variable_i,7]!=0){
                                f.con[f.con_i,A_variable[A_variable_i,7]]=-1;
			}

			f.con[f.con_i,length(local_l)+nrow(local_A)+SOS_i]=-1;
			SOS_i=SOS_i+1;
			f.rhs_i=f.rhs_i+1;
			f.con_i=f.con_i+1;

                        f.con[f.con_i,A_variable[A_variable_i,5]]=+1;
                        f.con[f.con_i,length(local_l)+nrow(local_A)+SOS_i]=-1;
			f.rhs[f.rhs_i]=0;
                        SOS_i=SOS_i+1;
                        f.rhs_i=f.rhs_i+1;
                        f.con_i=f.con_i+1;

                }else{
                        f.con[f.con_i,A_variable[A_variable_i,5]]=+1;
			f.rhs[f.rhs_i]=max(local_A[A_variable_i,4]-local_A[A_variable_i,6], 0);
                        f.rhs_i=f.rhs_i+1;
                        f.con_i=f.con_i+1;
                }
        }
        if(A_variable[A_variable_i,7]!=0){
                if(!(A_variable[A_variable_i,4]==0 && A_variable[A_variable_i,6]==0 && A_variable[A_variable_i,5]==0)){
                        f.con[f.con_i,A_variable[A_variable_i,7]]=+1;
                        f.rhs[f.rhs_i]=0;
                        if(A_variable[A_variable_i,4]==0){
                                f.rhs[f.rhs_i]=f.rhs[f.rhs_i]-local_A[A_variable_i,4];
                        }else{
                                f.con[f.con_i,A_variable[A_variable_i,4]]=+1;
                        }
                        if(A_variable[A_variable_i,6]==0){
                                f.rhs[f.rhs_i]=f.rhs[f.rhs_i]+local_A[A_variable_i,6];
                        }else{
                                f.con[f.con_i,A_variable[A_variable_i,6]]=-1;
                        }
                        if(A_variable[A_variable_i,5]!=0){
                                f.con[f.con_i,A_variable[A_variable_i,5]]=-1;
                        }

                        f.con[f.con_i,length(local_l)+nrow(local_A)+SOS_i]=-1;
                        SOS_i=SOS_i+1;
                        f.rhs_i=f.rhs_i+1;
                        f.con_i=f.con_i+1;

                        f.con[f.con_i,A_variable[A_variable_i,7]]=+1;
                        f.con[f.con_i,length(local_l)+nrow(local_A)+SOS_i]=-1;
                        f.rhs[f.rhs_i]=0;
			SOS_i=SOS_i+1;
                        f.rhs_i=f.rhs_i+1;
                        f.con_i=f.con_i+1;
                }else{
                        f.con[f.con_i,A_variable[A_variable_i,7]]=+1;
                        f.rhs[f.rhs_i]=max(local_A[A_variable_i,6]-local_A[A_variable_i,4], 0);
                        f.rhs_i=f.rhs_i+1;
                        f.con_i=f.con_i+1;
                }
        }
}

##################################
        if(nrow_A_CN!=0){
                for(l_i in 1:nrow_A_CN){
##                	if(!is.na(l[[l_i]][[5]][4]) && !l[[l_i]][[5]][4]>=ABSOLUTE_CN_limit ){
                        if(!is.na(l[[l_i]][[5]][4]) && !l[[l_i]][[5]][4]>=Inf ){
#				if(l[[l_i]][[5]][5]<=0){
#					which_modal_cn_delta=modal_cn_delta_for_zero_marker;
#				}else if(l[[l_i]][[5]][4]>ABSOLUTE_CN_limit){
#					which_modal_cn_delta=modal_cn_delta_for_high_copy+round(l[[l_i]][[5]][4]/10);
#				}else{
#					which_modal_cn_delta=modal_cn_delta;
#				}
                                if(l[[l_i]][[5]][4]>static_CN_limit){
                                        which_modal_cn_delta1=-modal_cn_delta_for_high_copy +round(l[[l_i]][[5]][4]/10);
                                        which_modal_cn_delta2=modal_cn_delta_for_high_copy +round(l[[l_i]][[5]][4]/10);
###############################################################################################
#				}else if(l[[l_i]][[5]][6]<seg.q.tab.threshold){
#################################### #############################################################
                                }else if(l[[l_i]][[5]][7]<seg.q.tab.threshold){

					if(l[[l_i]][[5]][4]<l[[l_i]][[5]][6]){
	                                        which_modal_cn_delta1=0;
	                                        which_modal_cn_delta2=(l[[l_i]][[5]][6]-l[[l_i]][[5]][4]);
					}else{
                                                which_modal_cn_delta1=(l[[l_i]][[5]][6]-l[[l_i]][[5]][4]);
                                                which_modal_cn_delta2=0;
					}
				}else{
                                        which_modal_cn_delta1=-modal_cn_delta_for_zero_marker;
                                        which_modal_cn_delta2=modal_cn_delta_for_zero_marker;
				}

				f.con[f.con_i, l_i]=1;
				f.rhs[f.rhs_i]=l[[l_i]][[5]][4]+which_modal_cn_delta1;
	                        f.rhs_i=f.rhs_i+1;
	                        f.con_i=f.con_i+1;
                              	f.con[f.con_i, l_i]=1;
                                f.rhs[f.rhs_i]=l[[l_i]][[5]][4]+which_modal_cn_delta2;
                                f.rhs_i=f.rhs_i+1;
                                f.con_i=f.con_i+1;
                        }
		}
                for(l_i in 1:nrow_A_CN){
##			if(!is.na(l[[l_i]][[5]][4]) && l[[l_i]][[5]][4]>=ABSOLUTE_CN_limit){
                        if(!is.na(l[[l_i]][[5]][4]) && l[[l_i]][[5]][4]>=Inf){
                                f.con[f.con_i, l_i]=1;
                                f.rhs[f.rhs_i]=l[[l_i]][[5]][4];
                                f.rhs_i=f.rhs_i+1;
                                f.con_i=f.con_i+1;
	                }
		}
        }

#################### for maximum of SVs 
	 SV_min=0;
         f.con.add=rep(0,ncol(f.con));
##         f.con=rbind(f.con,f.con.add);
         f.rhs[f.rhs_i]=SV_min;
         ##f.dir=c(f.dir, ">");
	 for(l_i in (nrow_A_CN+1):length(local_l)){
		f.con[f.con_i,l_i]=1;
  	 }

##	f.rhs_i=f.rhs_i+1;
##	f.con_i=f.con_i+1;
####################


#	for(l_i in 1:nrow_A_CN){
#		f.con[f.conf_i,l_i]=1;
#		f.rhs[f.rhs_i]=max(l[[l_i]][[1]]);
#		f.conf_i=f.conf_i+1;
#		f.rhs_i=f.rhs_i+1;
#	}
#	for(l_i in (nrow_A_CN+1):length(l)){
#		f.conf[f.conf_i,l_i]=1;
#		f.rhs[f.rhs_i]=min(l[[l_i]][[1]],);
#
#	}
	################################################################################## lpSolveApi
	bglp=list();
	bglp_i=1;
	bglp[[bglp_i]]=make.lp(nrow(f.con),length(f.obj))
	for(bglp_c in 1:ncol(f.con)){
		set.column(bglp[[bglp_i]], bglp_c, f.con[,bglp_c])
	}
        set.constr.type(bglp[[bglp_i]], f.dir)
	set.rhs(bglp[[bglp_i]], f.rhs)
	set.type(bglp[[bglp_i]], 1:length(f.obj), "integer");

	if(SV_SOS_count!=0){
		SOS_i=1;
		while(length(local_l)+nrow(local_A)+SOS_i+1 <=(length(local_l)+nrow(local_A)+SV_SOS_count)){
                        tryCatch({
				add.SOS(bglp[[bglp_i]],SOS_i,type=1,(SOS_i%/%2)+1,c(length(local_l)+nrow(local_A)+SOS_i,length(local_l)+nrow(local_A)+SOS_i+1),c(1,2))
	                        }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
			SOS_i=SOS_i+2;
		}
	}

        set.objfn(bglp[[bglp_i]], f.obj)############################## 반드시 맨 나중에 입력? $$$$$$$$$$$$$$$$

	lp.control(bglp[[bglp_i]],timeout=time_out)
	lpSolve_status=solve(bglp[[bglp_i]])
	#lp_objval=get.objective(bglp);
	lpSolve_out=get.variables(bglp[[bglp_i]]);
        lp_objval=sum(f.obj*lpSolve_out);
        print("lpSolve_status");
        print(lpSolve_status);
        print("lpSolve_objval");
        print(lp_objval);

	if(lpSolve_status!=0)
		print("lpSolve_error");
	if(lpSolve_status==7)
		print("timeout");

	prev_solution=lpSolve_out;
	#prev_objval=Inf;
	prev_objval=lp_objval;
	prev_status=lpSolve_status;
##        while(lp_objval <= prev_objval && !all(lpSolve_out == prev_solution)){
         while(lp_objval- prev_objval < lp_obj_thres && SV_min < max_SV_min &&( lpSolve_status == 0 || lpSolve_status == 1) && lp_iter_thres > 0){
##         while(lp_objval <= prev_objval && lpSolve_status==0){
        ## while(abs(lp_objval-prev_objval)<0.1 && lpSolve_status==0){
                prev_objval=lp_objval;
                prev_solution=lpSolve_out;
		prev_status=lpSolve_status;
##                SV_min=SV_min+sum(prev_solution[(nrow_A_CN+1):length(local_l)])+1;
                SV_min=sum(prev_solution[(nrow_A_CN+1):length(local_l)])+1;

		f.rhs[f.rhs_i]=SV_min;
		bglp_i=bglp_i+1;
	        bglp[[bglp_i]]=make.lp(nrow(f.con),length(f.obj))
        	for(bglp_c in 1:ncol(f.con)){
                	set.column(bglp[[bglp_i]], bglp_c, f.con[,bglp_c])
	       	}	
        	set.constr.type(bglp[[bglp_i]], f.dir)
	        set.rhs(bglp[[bglp_i]], f.rhs)
	        set.type(bglp[[bglp_i]], 1:length(f.obj), "integer");
		if(SV_SOS_count!=0){
                	SOS_i=1;
	                while(length(local_l)+nrow(local_A)+SOS_i+1 <=(length(local_l)+nrow(local_A)+SV_SOS_count)){
                        tryCatch({
	                        add.SOS(bglp[[bglp_i]],SOS_i,type=1,(SOS_i%/%2)+1,c(length(local_l)+nrow(local_A)+SOS_i,length(local_l)+nrow(local_A)+SOS_i+1),c(1,2))
                                }, error=function(e){cat("ERROR :", conditionMessage(e), "\n")})
                                SOS_i=SOS_i+2;
	                }
	        }

       		set.objfn(bglp[[bglp_i]], f.obj)############################## 반드시 맨 나중에 입력? $$$$$$$$$$$$$$$$
		print("SV_min");
		print(SV_min);
	        lp.control(bglp[[bglp_i]],timeout=time_out)
		lpSolve_status=solve(bglp[[bglp_i]])
#       		lp_objval=get.objective(bglp);
#		print(lp_objval)
	        lpSolve_out=get.variables(bglp[[bglp_i]]);
	        
		lp_objval=sum(f.obj*lpSolve_out);
		lp_iter_thres = lp_iter_thres - 1;

		print("lpSolve_status");
		print(lpSolve_status);
		print("lpSolve_objval");
		print(lp_objval);
	}
	lp_objval=prev_objval;
	lpSolve_out=prev_solution;
	print("FINAL lpSolve_status");
	print(prev_status);
	print("FINAL SV constraint satisfied");
	print(SV_constraint_explicit(local_A, A_variable, lpSolve_out));

#       lpSolve_out= lp_out$solution;
#	lp_objval=lp_out$objval;
	for(lpSolve_out_i in 1:length(local_l)){
		local_l[[lpSolve_out_i]][[3]]=lpSolve_out[lpSolve_out_i];
	}
        for(l_i in 1:length(local_l)){
                for(l_i_2 in 1:ncol(local_l[[l_i]][[2]])){
                        local_A[local_l[[l_i]][[2]][1,l_i_2],local_l[[l_i]][[2]][2,l_i_2]]=local_l[[l_i]][[3]];
                }
        }

        return(list(local_l,local_A,lp_objval));

}
Integer_prog_H_negative <- function(L, nrow_A_CN){
        local_l=L[[1]];
        local_A=L[[2]];
        f.obj=c();
#       f.con=matrix(0, nrow=nrow(A)+length(l), ncol=length(l));
        f.con=matrix(0, nrow=nrow(local_A), ncol=length(local_l));
        f.con_i=1;
        f.dir=rep("<", nrow(f.con));
        f.rhs=c();
        f.rhs_i=1;
        A_variable=local_A;
        A_variable[1:nrow(local_A),4:7]=0;
        for(l_i in 1:length(local_l)){
                f.obj[l_i]=0;
                for(l_i_2 in 1:ncol(local_l[[l_i]][[2]])){
                        if(local_l[[l_i]][[2]][2,l_i_2]== 4 || local_l[[l_i]][[2]][2,l_i_2]== 7){
                                f.obj[l_i]=f.obj[l_i]-1;
                        }else{
                                f.obj[l_i]=f.obj[l_i]+1;
                        }
                        A_variable[local_l[[l_i]][[2]][1,l_i_2],local_l[[l_i]][[2]][2,l_i_2]]=l_i;
                }
        }
        for(A_i in 1:nrow(A_variable)){
                f.rhs[f.rhs_i]=0;
                if(A_variable[A_i,4]==0){
                        f.rhs[f.rhs_i]=f.rhs[f.rhs_i]-local_A[A_i,4];
                }else{
                        f.con[f.con_i,A_variable[A_i,4]]=f.con[f.con_i,A_variable[A_i,4]]+1;
                }
                if(A_variable[A_i,5]==0){
                        f.rhs[f.rhs_i]=f.rhs[f.rhs_i]+local_A[A_i,5];
                }else{
                        f.con[f.con_i,A_variable[A_i,5]]=f.con[f.con_i,A_variable[A_i,5]]-1;
                }
                if(A_variable[A_i,6]==0){
                        f.rhs[f.rhs_i]=f.rhs[f.rhs_i]+local_A[A_i,6];
                }else{
                        f.con[f.con_i,A_variable[A_i,6]]=f.con[f.con_i,A_variable[A_i,6]]-1;
                }
                if(A_variable[A_i,7]==0){
                        f.rhs[f.rhs_i]=f.rhs[f.rhs_i]-local_A[A_i,7];
                }else{
                        f.con[f.con_i,A_variable[A_i,7]]=f.con[f.con_i,A_variable[A_i,7]]+1;
                }
                f.rhs_i=f.rhs_i+1;
                f.con_i=f.con_i+1;
        }        
	
	for(A_i in 1:nrow(A_variable)){
                if(A_variable [A_i,5] !=0){
                        f.con.add = rep(0, length(local_l));
                        f.dir=c(f.dir, ">");
                        if(A_variable[A_i,4] == 0){
                                f.rhs[f.rhs_i]=-local_A[A_i,4];
                                f.con.add[A_variable [A_i,5]]=-1;
                        }else{
                                f.rhs[f.rhs_i]=0;
                                f.con.add[A_variable [A_i,4]]=1;
                                f.con.add[A_variable [A_i,5]]=-1;
                        }
                        f.rhs_i=f.rhs_i+1;
                        f.con=rbind(f.con,f.con.add);
                }
                if(A_variable [A_i,7] !=0){
                        f.con.add = rep(0, length(local_l));
                        f.dir=c(f.dir, ">");
                        if(A_variable[A_i,6] == 0){
                                f.rhs[f.rhs_i]=-local_A[A_i,6];
                                f.con.add[A_variable [A_i,7]]=-1;
                        }else{
                                f.rhs[f.rhs_i]=0;
                                f.con.add[A_variable [A_i,6]]=1;
                                f.con.add[A_variable [A_i,7]]=-1;
                        }
                        f.rhs_i=f.rhs_i+1;
                        f.con=rbind(f.con,f.con.add);
                }
        }

#       for(l_i in 1:nrow_A_CN){
#               f.con[f.conf_i,l_i]=1;
#               f.rhs[f.rhs_i]=max(l[[l_i]][[1]]);
#               f.conf_i=f.conf_i+1;
#               f.rhs_i=f.rhs_i+1;
#       }
#       for(l_i in (nrow_A_CN+1):length(l)){
#               f.conf[f.conf_i,l_i]=1;
#               f.rhs[f.rhs_i]=min(l[[l_i]][[1]],);
#
#       }
        lp_out=lp("min", f.obj, f.con, f.dir, f.rhs);
        lpSolve_out= lp_out$solution;
        lp_objval=lp_out$objval;

        for(lpSolve_out_i in 1:length(lpSolve_out)){
                local_l[[lpSolve_out_i]][[3]]=lpSolve_out[lpSolve_out_i];
        }
        for(l_i in 1:length(local_l)){
                for(l_i_2 in 1:ncol(local_l[[l_i]][[2]])){
                        local_A[local_l[[l_i]][[2]][1,l_i_2],local_l[[l_i]][[2]][2,l_i_2]]=local_l[[l_i]][[3]];
                }
        }

        return(list(local_l,local_A,lp_objval));

}




Integer_prog<- function (list, l_index,l_max, nrow_A_CN){
	l=list[[1]];
	A=list[[2]];
        prev_nochange_A=A;
        prev_nochange_l=l;

	if(l_index<=nrow_A_CN){
##		prev_nochange_A=A;
##		prev_nochange_l=l;
		for(i in 1:length(l[[l_index]][[1]])){
		        prev_A=A;
			prev_l=l;
		        prev_obj=Obj_func(prev_A);
			A=prev_nochange_A;
			l=prev_nochange_l;

                        for(target_number in 1:ncol(l[[l_index]][[2]])){
                                A[l[[l_index]][[2]][1,target_number], l[[l_index]][[2]][2,target_number]]=l[[l_index]][[1]][i];
                        }

			current_CN=l[[l_index]][[1]][i];
                        l[[l_index]][[3]]=current_CN;
			if(l_index<nrow_A_CN&&l[[l_index]][[5]][1]==l[[l_index+1]][[5]][1]&&l[[l_index]][[5]][3]+1==l[[l_index+1]][[5]][2]){
				CN_intersection=intersect(l[[l_index+1]][[1]],c(current_CN-1, current_CN, current_CN+1));	
				if(length(CN_intersection)!=0){
					l[[l_index+1]][[1]]=CN_intersection;   ##adjacent CN constraint##
				}else{
					A[1,4]=10000 ### The case is absolutely impossible  10000~infinite number##

				}
			}
			##print(l);
			if(A[1,4]!=10000){
				Integer_prog_return=Integer_prog(list(l,A), l_index+1,l_max, nrow_A_CN);	
				l=Integer_prog_return[[1]];
				A=Integer_prog_return[[2]];
			}
			##print(A);
			##print(Obj_func(A));
			if(Obj_func(A)>=prev_obj){
				A=prev_A;
				l=prev_l;
			}
				
		}
	}
	if(l_index>nrow_A_CN && l_index <= l_max){
		if(l[[l_index]][[4]]==1){
			remain_degree1=A[l[[l_index]][[2]][1,1], l[[l_index]][[2]][2,1]-1]-A[l[[l_index]][[2]][1,1], 12-l[[l_index]][[2]][2,1]-1];
			remain_degree2=A[l[[l_index]][[2]][1,2], l[[l_index]][[2]][2,2]-1]-A[l[[l_index]][[2]][1,2], 12-l[[l_index]][[2]][2,2]-1];
		 	SV_number_decided = max(0,min(remain_degree1, remain_degree2));
                        for(target_number in 1:ncol(l[[l_index]][[2]])){
                                A[l[[l_index]][[2]][1,target_number], l[[l_index]][[2]][2,target_number]]=SV_number_decided;
                        }
			l[[l_index]][[3]]=SV_number_decided;
                        Integer_prog_return=Integer_prog(list(l,A), l_index+1,l_max, nrow_A_CN);
                        l=Integer_prog_return[[1]];
                        A=Integer_prog_return[[2]];

		}else{
                for(i in 1:length(l[[l_index]][[1]])){
                        prev_A=A;
                        prev_l=l;
                        prev_obj=Obj_func(prev_A);
			A=prev_nochange_A;
                        l=prev_nochange_l;


                        for(target_number in 1:ncol(l[[l_index]][[2]])){
                                A[l[[l_index]][[2]][1,target_number], l[[l_index]][[2]][2,target_number]]=l[[l_index]][[1]][i];
                        }

                        current_CN=l[[l_index]][[1]][i];
                        l[[l_index]][[3]]=current_CN;


                        Integer_prog_return=Integer_prog(list(l,A), l_index+1,l_max, nrow_A_CN);
                        l=Integer_prog_return[[1]];
                        A=Integer_prog_return[[2]];
                        if(Obj_func(A)>=prev_obj){
                                A=prev_A;
                                l=prev_l;
                        }

                }


		}
	}


	return(list(l,A));
}

Obj_func<- function (A){
sum=0;
for(i in 1:nrow(A)){
	sum=sum+abs(A[i,6]-A[i,7]-(A[i,4]-A[i,5]));
}
return(sum);
}

while(nrow(SV_set)!=0){
##for(i in 7:7){
	i=1;
	SV=SV_set;
	A=matrix(, ncol=7, nrow=0);
#	A_CN=matrix(, ncol=4, nrow=0);
        A_CN=matrix(, ncol=7, nrow=0);

	SV_line=SV[i,];
#	SV=SV[-i,];
#	related_set=related_set_find(SV_line,SV,A);
	related_set=related_set_find_BFS_ver(SV_line,SV,A);
	A=related_set[[1]];
	print(A);
	A_SV_remain=related_set[[2]];
	
	for(j in 1:nrow(A)){
		if(is.na(new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],5])!=T&&new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],4]>marker_size&& new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],"expected_cn"]!=Inf&&!(abs(new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],"expected_cn"]-new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],"modal_cn"])>min_expected_modal_cn_diff)){
			A[j,4]=new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],18]}else{
			A_CN=rbind(A_CN,c(new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],1], new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],2], new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],3]	, new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],18],new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],"Probes"],new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],"alternative_modal_cn"],new[new$Chromosome==A[j,1]&new$End.bp==A[j,2],"seg.q.tab"]));
			}

		if(is.na(new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],5])!=T&&new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],4]>marker_size&& new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],"expected_cn"]!=Inf &&!(abs(new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],"expected_cn"]-new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],"modal_cn"])>min_expected_modal_cn_diff)){
			A[j,6]=new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],18]}else{
                        A_CN=rbind(A_CN,c(new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],1], new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],2], new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],3]  , new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],18],new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],"Probes"],new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],"alternative_modal_cn"],new[new$Chromosome==A[j,1]&new$Start.bp==A[j,3],"seg.q.tab"]));
			}
	}
	A_CN=unique(A_CN);
	A=A[order(A[,1],A[,2]),]
	
	if(nrow(A_CN)>1){
		A_CN=A_CN[order(A_CN[,1],A_CN[,2]),];
	}
	##print(SV_set);
##	print(A_SV_remain);
#	A_SV=SV_matrix_minus(SV_set,A_SV_remain);
        A_SV=SV_set[related_set[[3]],];

	SV_set=A_SV_remain;

	l<-list();
	l_index=1;

        if(nrow(A_CN)!=0){
        for(CN_index in 1:nrow(A_CN)){

##### SV target find #####
		SV_target=c();
		for(A_SV_i in 1:nrow(A_SV)){
                        if(A_SV[A_SV_i,2]==A_CN[CN_index,1]&&A_SV[A_SV_i,3]==A_CN[CN_index,2]-1)
                                SV_target=c(SV_target,A_SV_i);
                        if(A_SV[A_SV_i,4]==A_CN[CN_index,1]&&A_SV[A_SV_i,5]==A_CN[CN_index,2]-1)
                                SV_target=c(SV_target,A_SV_i);

			if(A_SV[A_SV_i,2]==A_CN[CN_index,1]&&A_SV[A_SV_i,3]==A_CN[CN_index,2])
				SV_target=c(SV_target,A_SV_i);
                        if(A_SV[A_SV_i,4]==A_CN[CN_index,1]&&A_SV[A_SV_i,5]==A_CN[CN_index,2])
                                SV_target=c(SV_target,A_SV_i);
		}		
############################

                target=matrix(,ncol=0,nrow=2)
                for(A_index in 1:nrow(A)){
                        if(A_CN[CN_index,1]==A[A_index,1]&&A_CN[CN_index,2]==A[A_index,3])
                                target=cbind(target,c(A_index,6));
                        if(A_CN[CN_index,1]==A[A_index,1]&&A_CN[CN_index,3]==A[A_index,2])
                                target=cbind(target,c(A_index,4));
                }
                if(is.na(A_CN[CN_index,4])==T){
			if(nrow(new[new$Chromosome==A_CN[CN_index,1]&new$End.bp==A_CN[CN_index,2]-1,])!=0){
				CN_l_modal=new[new$Chromosome==A_CN[CN_index,1]&new$End.bp==A_CN[CN_index,2]-1,18];
				m_l_size=new[new$Chromosome==A_CN[CN_index,1]&new$End.bp==A_CN[CN_index,2]-1,4];
				if(is.na(CN_l_modal)!=T&&m_l_size>marker_size){
					CN_l_range=c(CN_l_modal, CN_l_modal-1, CN_l_modal+1)
				}else{
					CN_l_range=c(0,1,2,3,4,5,6,7,8,9,10);
				}
			}else{
				CN_l_range=c(0,1,2,3,4,5,6,7,8,9,10);
			}

			if(nrow(new[new$Chromosome==A_CN[CN_index,1]&new$Start.bp==A_CN[CN_index,3]+1,])!=0){
	                        CN_r_modal=new[new$Chromosome==A_CN[CN_index,1]&new$Start.bp==A_CN[CN_index,3]+1,18];
				m_r_size=new[new$Chromosome==A_CN[CN_index,1]&new$Start.bp==A_CN[CN_index,3]+1,4];
	                        if(is.na(CN_r_modal)!=T&&m_r_size>marker_size){
	                                CN_r_range=c(CN_r_modal, CN_r_modal-1, CN_r_modal+1)
	                        }else{
					CN_r_range=c(0,1,2,3,4,5,6,7,8,9,10);
				}
			}else{
				CN_r_range=c(0,1,2,3,4,5,6,7,8,9,10);
			}

			if(length(intersect(CN_l_range, CN_r_range))!=0){
                       		l[[l_index]]<-list(intersect(CN_l_range, CN_r_range),target,0, SV_target+nrow(A_CN), A_CN[CN_index,] )
			}else{
				l[[l_index]]<-list(c(0,1,2,3,4,5,6,7,8,9,10),target,0, SV_target+nrow(A_CN), A_CN[CN_index,] )
			}
                }else{
                        l[[l_index]]<-list(c(A_CN[CN_index,4], A_CN[CN_index,4]-1, A_CN[CN_index,4]+1),target,0, SV_target+nrow(A_CN), A_CN[CN_index,])
                }



                l_index=l_index+1;

        }
        }

       for(SV_index in 1:nrow(A_SV)){
                target=matrix(,ncol=0,nrow=2)

                for(A_index in 1:nrow(A)){
                        if(A_SV[SV_index,2]==A[A_index,1]&&A_SV[SV_index,3]==A[A_index,2])
                                target=cbind(target,c(A_index,5));
                        if(A_SV[SV_index,2]==A[A_index,1]&&A_SV[SV_index,3]==A[A_index,3])
                                target=cbind(target,c(A_index,7));
                        if(A_SV[SV_index,4]==A[A_index,1]&&A_SV[SV_index,5]==A[A_index,2])
                                target=cbind(target,c(A_index,5));
                        if(A_SV[SV_index,4]==A[A_index,1]&&A_SV[SV_index,5]==A[A_index,3])
                                target=cbind(target,c(A_index,7));

                }




	        l[[l_index]]<-list(c(0,1,2),target,0);
                l_index=l_index+1;
        }

	l_SV_start=nrow(A_CN)+1;
	for(l_SV_index in l_SV_start:length(l)){
		count =0;
		for( ll in l_SV_start:length(l))
			if(l[[l_SV_index]][[2]][1,1]==l[[ll]][[2]][1,1]||l[[l_SV_index]][[2]][1,1]==l[[ll]][[2]][1,2]||l[[l_SV_index]][[2]][1,2]==l[[ll]][[2]][1,1]||l[[l_SV_index]][[2]][1,2]==l[[ll]][[2]][1,2])
			count=count+1;


		l[[l_SV_index]][[4]]=count;  ## count whether the other SV exists in two breakpoints of the SV.
		print("count");
		print(count);

  	}
##	print(l);
	if(length(l)<1){
		Integer_prog_return=Integer_prog(list(l,A),1,length(l), nrow(A_CN));
        	l=Integer_prog_return[[1]];
	        A=Integer_prog_return[[2]];

	}else{
		Integer_prog_H_positive_return=Integer_prog_H_positive(list(l,A),nrow(A_CN));	
              #  Integer_prog_H_negative_return=Integer_prog_H_negative(list(l,A),nrow(A_CN));
	#	if(Integer_prog_H_positive_return[[3]]<Integer_prog_H_negative_return[[3]]){
	        	l=Integer_prog_H_positive_return[[1]];
		        A=Integer_prog_H_positive_return[[2]];
	#	}else{
         #               l=Integer_prog_H_negative_return[[1]];
          #              A=Integer_prog_H_negative_return[[2]];
	#	}

	}

##	l=Integer_prog_return[[1]];
##	A=Integer_prog_return[[2]];
	A_SV=cbind(A_SV, CN=-1);

        if(nrow(A_CN)!=0){
        for(A_CN_i in 1:nrow(A_CN)){
                A_CN[A_CN_i,4]=l[[A_CN_i]][[3]];
                new[new$Chromosome==A_CN[A_CN_i,1]&new$Start.bp==A_CN[A_CN_i,2]&new$End.bp==A_CN[A_CN_i,3],18]=A_CN[A_CN_i,4];
        }
        }else{
		A_CN_i=0
	};

	
################# false positive SV test ########################
	false_SVs=c();
	for(A_SV_i in 1:nrow(A_SV)){
		print(A_SV[A_SV_i,]);
		false_SV_test=0;
		if(strsplit(as.character(A_SV[A_SV_i,6]),"to")[[1]][1]=="5"){where1="Start.bp";where2_d=-1;}else{where1="End.bp";where2_d=+1;};
		where1=which(new[,where1]==A_SV[A_SV_i,3] & new[,1]==A_SV[A_SV_i,2]);
		where2=where1+where2_d;
		if(is.na(new[where1,"modal_cn"]) ||  is.na(new[where2,"modal_cn"]) ||is.na(new[where1,"n_probes"]) ||is.na(new[where2,"n_probes"]) ||  new[where1,"n_probes"] <false_marker_size || new[where2,"n_probes"] < false_marker_size){
			if(is.na(new[where1,"raw_expected_cn"]) || is.na(new[where2,"raw_expected_cn"]) || max(0,round(new[where1,"raw_expected_cn"])) <= max(0,round(new[where2,"raw_expected_cn"]))){
				false_SV_test=false_SV_test+1;
			}
		}
                if(strsplit(as.character(A_SV[A_SV_i,6]),"to")[[1]][2]=="5"){where1="Start.bp";where2_d=-1;}else{where1="End.bp";where2_d=+1;};
                where1=which(new[,where1]==A_SV[A_SV_i,5] & new[,1]==A_SV[A_SV_i,4]);
                where2=where1+where2_d;
                if(is.na(new[where1,"modal_cn"]) ||  is.na(new[where2,"modal_cn"]) ||is.na(new[where1,"n_probes"]) ||is.na(new[where2,"n_probes"]) ||  new[where1,"n_probes"] <false_marker_size || new[where2,"n_probes"] < false_marker_size){
                        if(is.na(new[where1,"raw_expected_cn"]) || is.na(new[where2,"raw_expected_cn"]) || max(0,round(new[where1,"raw_expected_cn"])) <= max(0,round(new[where2,"raw_expected_cn"]))){
	                        false_SV_test=false_SV_test+1;
			}
		}
		
		if(false_SV_test!=0){
			false_SVs=c(false_SVs, 1);
		}else{
			false_SVs=c(false_SVs, 0);
		}
	}

#############################################
	for(A_SV_i in 1:nrow(A_SV)){
		if(false_SVs[A_SV_i]==1){
       ##         if(0){
			A_SV[A_SV_i,15]=0;
		}else{
                        A_SV[A_SV_i,15]=l[[A_CN_i+A_SV_i]][[3]];
		}
		SV_table[SV_table$V2==A_SV[A_SV_i,2]&SV_table$V3==A_SV[A_SV_i,3]&SV_table$V4==A_SV[A_SV_i,4]&SV_table$V5==A_SV[A_SV_i,5],15]=A_SV[A_SV_i,15];
	}

	print(A);
	print(A_SV);
	print(A_CN);
	print("Obj");
	print(Obj_func(A));
}

SV_table[which(SV_table[,2]==23),2]<-"X"
SV_table[which(SV_table[,4]==23),4]<-"X"
##levels(SV_table[,2])[levels(SV_table[,2])=="23"]<-"X"
##levels(SV_table[,4])[levels(SV_table[,4])=="23"]<-"X"

SV_table_filtered=SV_table[SV_table$V15!=0,];
SV_table_filtered=SV_table_filtered[,-15];

write.table(SV_table, paste(args[3],".CN_opt",sep=""), quote=F, row.names=F, col.names=F, sep="\t")
write.table(SV_table_filtered, paste(args[3],".CN_opt.filtered",sep=""), quote=F, row.names=F, col.names=F,sep="\t")

options(scipen = 999)
write.table(new, "copy_numbers.CN_opt", quote=F, row.names=F, sep="\t")

############################

for(i in 1:nrow(new)){
	if(length(which(ABS_output_duplicated$Chromosome==new[i,"Chromosome"]&ABS_output_duplicated$Start.bp==new[i,"Start.bp"] & ABS_output_duplicated$End.bp==new[i,"End.bp"]))==1){
		new[i,"expected_cn"]=ABS_output_duplicated[which(ABS_output_duplicated$Chromosome==new[i,"Chromosome"]&ABS_output_duplicated$Start.bp==new[i,"Start.bp"] & ABS_output_duplicated$End.bp==new[i,"End.bp"]),"expected_cn"];
	}
}
write.table(new, "copy_numbers.CN_opt.debug", quote=F, row.names=F, sep="\t")


################################################
