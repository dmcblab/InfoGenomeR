args=commandArgs(T)
min_dm_copy=5

library(plyr)
path=scan("euler_paths.0", what="character", sep="\n",skip=1)
for( i in 1:length(path)){
	path[i] = strsplit(as.character(path[i]), "\t")
}


p_eq= function(lp1, lp2){
	if(length(lp1) == length(lp2)){
		check1=T;
		check2=T;
		for(lp1_i in 1:length(lp1)){
			if(lp1[[lp1_i]] != lp2[[lp1_i]]){
				check1=F
			}
			if(lp1[[lp1_i]] != lp2[[length(lp2)-lp1_i+1]]){
				check2=F
			}
		}
		if(any(check1, check2)){
			return(T)
		}else{
			return(F)
		}

	}else{
		return(F)
	}
}


entropy= function(lp){
        lp_c=c()
        for(lpi in 1:length(lp)){
                if(lpi!=1){
                        for(lpj in 1:(lpi-1)){
                                if(p_eq(lp[[lpj]] , lp[[lpi]])){
                                        lp_c[lpi]=lp_c[lpj]
                                }
                        }
                        if(length(lp_c)<lpi){
                                        lp_c[lpi]=max(lp_c)+1
                        }
                }else{
                        lp_c=c(1)
                }
        }
	w=(lp!="-2")
	lp_c=lp_c[w]
        lp_ent=count(lp_c)
        return(sum(-(lp_ent[,2]/sum(lp_ent[,2]))*log(lp_ent[,2]/sum(lp_ent[,2]),base=2)))

}



op=path

path_count=rep(0, length(path))
for(i in 1:length(path)){
	e_path_count=0
	for(j in 1:length(path)){
		if(length(path[[i]]) == length(path[[j]]) && (all(path[[i]] == path [[j]]) || all (rev(path[[i]]) == path[[j]]))){
			e_path_count=e_path_count+1
		}
	}
	path_count[i]=e_path_count;
}

#v=c()
for( i in 1:length(op)){
	if (path[[i]][1] == "-1"){

		ent1= entropy(op)

		n=path[[i]][2]
		ns=path[[i]][3]
		ne=path[[i]][length(path[[i]])-2]
		

		ij_l=c()
		ik_l=c()

		for( j in 1:length(path)){
			for( k in 1:length(path[[j]])  ){
				if(path[[j]][1] == "-1"){
					next
				}
				if(path[[j]][k] == n){
					ij_l=c(ij_l,j)
					ik_l=c(ik_l,k)
				}
			}
		}

		if(length(ij_l)>0){
			ent2=c()

			alt_path=list()
			for(li in 1:length(ij_l)){

				path=op
				ij=ij_l[li]
				ik=ik_l[li]
				p=c()
				if(ik > 1 && ik < length(path[[ij]]) && path[[ij]][ik-1]  == ne){
					p=path[[i]][3:(length(path[[i]])-1)]
				}else if(ik> 1 && ik < length(path[[ij]]) && path[[ij]][[ik+1]] == ns){
					p=path[[i]][3:(length(path[[i]])-1)]
				}else if( ik > 1&& ik < length(path[[ij]]) && path[[ij]][ik-1] == ns){
					p=path[[i]][2:(length(path[[i]])-2)]
					p=rev(p)
				}else if(ik > 1 && ik < length(path[[ij]]) && path[[ij]][ik+1]==ne){
					p=path[[i]][2:(length(path[[i]])-2)]
					p=rev(p)
				}else{
					next;
				}
				if(length(p)>0){
					path[[ij]]=append(path[[ij]], p, after=ik)
					path[[i]]=-2;
				}
				ent2= c(ent2,entropy(path))
				alt_path[[length(alt_path)+1]]=path
			}
			if(args[1]!="dup" && path_count[i]>min_dm_copy){
				if(ent1 > min(ent2)){
					w=which.min(ent2)
					op=alt_path[[w]]
				}
			}else{
				if(length(ent2)>0){
					w=which.min(ent2)
					op=alt_path[[w]]
				}
			}
		}
	}
}


path=op
v=which(path=="-2")
path[v]=c()


for( i in 1:length(path)){


	if(path[[i]][1] == "-1"){
		if( abs(as.numeric(path[[i]][3]) - as.numeric(path[[i]][2])) ==1 ){
#		if( as.numeric(path[[i]][3]) == as.numeric(path[[i]][2]) -1){
			path[[i]]=c("circular",rep(path[[i]][2:(length(path[[i]])-2)],3))
		}else if(as.numeric(path[[i]][length(path[[i]])-1])-1 == as.numeric(path[[i]][length(path[[i]])-2])){
			 path[[i]]=c("circular",rep(path[[i]][3:(length(path[[i]])-1)],3))
		}
	}else{

		path[[i]] = c("linear", path[[i]])
	}
}

for( i in 1:length(path)){
	cat(paste(path[[i]],collapse="\t"))
	cat("\n")
}
