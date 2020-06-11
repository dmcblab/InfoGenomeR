AS_SV=read.table("SVs.CN_opt",stringsAsFactors=F)
AS_CN=read.table("copy_numbers.CN_opt.ACN.phased", header=T,stringsAsFactors=F)
AS_CN[AS_CN[,1]=="23",1]<-"X"

AS_SV=cbind(AS_SV, l_allele1=NA, l_allele2=NA, r_allele1=NA, r_allele2=NA, l_cn=NA, r_cn=NA);
new_AS_SV=AS_SV[0,];


for(i in 1:nrow(AS_SV)){
        f=strsplit(AS_SV[i,6],split="to");
	if(f[[1]][1]==3){
		la1=AS_CN[AS_CN$Chromosome==AS_SV[i,2]&AS_CN$End.bp==AS_SV[i,3],28];
		la2=AS_CN[AS_CN$Chromosome==AS_SV[i,2]&AS_CN$End.bp==AS_SV[i,3],29];
		la3=AS_CN[AS_CN$Chromosome==AS_SV[i,2]&AS_CN$Start.bp==AS_SV[i,3]+1,28];
		la4=AS_CN[AS_CN$Chromosome==AS_SV[i,2]&AS_CN$Start.bp==AS_SV[i,3]+1,29];
		lcn1=AS_CN[AS_CN$Chromosome==AS_SV[i,2]&AS_CN$End.bp==AS_SV[i,3],"modal_cn"];
		lcn2=AS_CN[AS_CN$Chromosome==AS_SV[i,2]&AS_CN$Start.bp==AS_SV[i,3]+1,"modal_cn"];
	}else{
		la1=AS_CN[AS_CN$Chromosome==AS_SV[i,2]&AS_CN$Start.bp==AS_SV[i,3],28];
		la2=AS_CN[AS_CN$Chromosome==AS_SV[i,2]&AS_CN$Start.bp==AS_SV[i,3],29];
                la3=AS_CN[AS_CN$Chromosome==AS_SV[i,2]&AS_CN$End.bp==AS_SV[i,3]-1,28];
		la4=AS_CN[AS_CN$Chromosome==AS_SV[i,2]&AS_CN$End.bp==AS_SV[i,3]-1,29];
		lcn1=AS_CN[AS_CN$Chromosome==AS_SV[i,2]&AS_CN$Start.bp==AS_SV[i,3],"modal_cn"];
		lcn2=AS_CN[AS_CN$Chromosome==AS_SV[i,2]&AS_CN$End.bp==AS_SV[i,3]-1,"modal_cn"];
	}

	 if(f[[1]][2]==3){
		ra1=AS_CN[AS_CN$Chromosome==AS_SV[i,4]&AS_CN$End.bp==AS_SV[i,5],28];
		ra2=AS_CN[AS_CN$Chromosome==AS_SV[i,4]&AS_CN$End.bp==AS_SV[i,5],29];
		ra3=AS_CN[AS_CN$Chromosome==AS_SV[i,4]&AS_CN$Start.bp==AS_SV[i,5]+1,28];
		ra4=AS_CN[AS_CN$Chromosome==AS_SV[i,4]&AS_CN$Start.bp==AS_SV[i,5]+1,29];
		rcn1=AS_CN[AS_CN$Chromosome==AS_SV[i,4]&AS_CN$End.bp==AS_SV[i,5],"modal_cn"];
		rcn2=AS_CN[AS_CN$Chromosome==AS_SV[i,4]&AS_CN$Start.bp==AS_SV[i,5]+1,"modal_cn"];
	}else{
		ra1=AS_CN[AS_CN$Chromosome==AS_SV[i,4]&AS_CN$Start.bp==AS_SV[i,5],28];
		ra2=AS_CN[AS_CN$Chromosome==AS_SV[i,4]&AS_CN$Start.bp==AS_SV[i,5],29];
		ra3=AS_CN[AS_CN$Chromosome==AS_SV[i,4]&AS_CN$End.bp==AS_SV[i,5]-1,28];
		ra4=AS_CN[AS_CN$Chromosome==AS_SV[i,4]&AS_CN$End.bp==AS_SV[i,5]-1,29];
		rcn1=AS_CN[AS_CN$Chromosome==AS_SV[i,4]&AS_CN$Start.bp==AS_SV[i,5],"modal_cn"];
		rcn2=AS_CN[AS_CN$Chromosome==AS_SV[i,4]&AS_CN$End.bp==AS_SV[i,5]-1,"modal_cn"];
	}

	if(is.na(la1) && is.na(ra1)){
	       new_AS_SV_i=nrow(new_AS_SV)+1;
		for(j in 1:15){
			new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
		}
		new_AS_SV[new_AS_SV_i,"l_cn"]=min(lcn1-lcn2  , rcn1-rcn2);
                new_AS_SV[new_AS_SV_i,"r_cn"]=min(lcn1-lcn2  , rcn1-rcn2);
	}else if(is.na(la1) && !is.na(ra1)){
		l_degrees=lcn1-lcn2;
		if(!is.na(ra3)){
			r_degree1=max(0,ra1-ra3);
			r_degree2=max(0,ra2-ra4);
			if(l_degrees <= max(r_degree1, r_degree2)){
			       new_AS_SV_i=nrow(new_AS_SV)+1;
				for(j in 1:15){
					new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
				}
				new_AS_SV[new_AS_SV_i,"l_cn"]=l_degrees;
				new_AS_SV[new_AS_SV_i,17+which.max(c(max(0,ra1-ra3), max(0,ra2-ra4)))]=l_degrees;
			}else{
                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i,"l_cn"]=max(r_degree1, r_degree2);
                                new_AS_SV[new_AS_SV_i,17+which.max(c(max(0,ra1-ra3), max(0,ra2-ra4)))]=max(r_degree1, r_degree2);
				l_degrees=l_degrees-max(r_degree1, r_degree2);

                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
				new_AS_SV[new_AS_SV_i, "l_cn"]=min(l_degrees, min(r_degree1, r_degree2));
				new_AS_SV[new_AS_SV_i, 17+which.min(rev(c(max(0,ra1-ra3), max(0,ra2-ra4))))]=min(l_degrees, min(r_degree1, r_degree2));

			}
		}else{
			if(l_degrees <=max(ra1,ra2)){
                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i,"l_cn"]=l_degrees;
                                new_AS_SV[new_AS_SV_i,17+which.max(c(ra1,ra2))]=l_degrees;
			}else{
                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i,"l_cn"]=max(ra1,ra2);
                                new_AS_SV[new_AS_SV_i,17+which.max(c(ra1,ra2))]=max(ra1,ra2);
				l_degrees=l_degrees-max(ra1,ra2);
                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i, "l_cn"]=min(l_degrees, min(ra1,ra2));
                                new_AS_SV[new_AS_SV_i, 17+which.min(rev(c(ra1,ra2)))]=min(l_degrees, min(ra1,ra2));
			}
		}
	}else if(!is.na(la1) && is.na(ra1)){
                r_degrees=rcn1-rcn2;
                if(!is.na(la3)){
                        l_degree1=max(0,la1-la3);
                        l_degree2=max(0,la2-la4);
                        if(r_degrees <= max(l_degree1, l_degree2)){
                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i,"r_cn"]=r_degrees;
                                new_AS_SV[new_AS_SV_i,15+which.max(c(max(0,la1-la3), max(0,la2-la4)))]=r_degrees;
                        }else{
                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i,"r_cn"]=max(l_degree1, l_degree2);
                                new_AS_SV[new_AS_SV_i,15+which.max(c(max(0,la1-la3), max(0,la2-la4)))]=max(l_degree1, l_degree2);
                                r_degrees=r_degrees-max(l_degree1, l_degree2);

                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i, "r_cn"]=min(r_degrees, min(l_degree1, l_degree2));
                                new_AS_SV[new_AS_SV_i, 15+which.min(rev(c(max(0,la1-la3), max(0,la2-la4))))]=min(r_degrees, min(l_degree1, l_degree2));

                        }
                }else{
                        if(r_degrees <=max(la1,la2)){
                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i,"r_cn"]=r_degrees;
                                new_AS_SV[new_AS_SV_i,15+which.max(c(la1,la2))]=r_degrees;
                        }else{
                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i,"r_cn"]=max(la1,la2);
                                new_AS_SV[new_AS_SV_i,15+which.max(c(la1,la2))]=max(la1,la2);
                                r_degrees=r_degrees-max(la1,la2);
                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i, "r_cn"]=min(r_degrees, min(la1,la2));
                                new_AS_SV[new_AS_SV_i, 15+which.min(rev(c(la1,la2)))]=min(r_degrees, min(la1,la2));
                        }
                }
	}else{
		if(!is.na(la3)){
			l_degree1=max(0,la1-la3);
			l_degree2=max(0,la2-la4);
		}else{
			l_degree=lcn1-lcn2;
			if(la1<=la2){
				if(l_degree <= la2){
					l_degree1=0;
					l_degree2=l_degree;
				}else{
					l_degree2=la2;
					l_degree=l_degree-la2;
					l_degree1=min(l_degree, la1);
				}
			}else{
                               if(l_degree <= la1){
                                        l_degree2=0;
                                        l_degree1=l_degree;
                                }else{
                                        l_degree1=la1;
                                        l_degree=l_degree-la1;
                                        l_degree2=min(l_degree, la2);
                                }
			}
		}

                if(!is.na(ra3)){
                        r_degree1=max(0,ra1-ra3);
                        r_degree2=max(0,ra2-ra4);
                }else{
                        r_degree=rcn1-rcn2;
                        if(ra1<=ra2){
                                if(r_degree <=ra2){
                                        r_degree1=0;
                                        r_degree2=r_degree;
                                }else{
                                        r_degree2=ra2;
                                        r_degree=r_degree-ra2;
                                        r_degree1=min(r_degree, ra1);
                                }
                        }else{
                               if(r_degree <= ra1){
                                        r_degree2=0;
                                        r_degree1=r_degree;
                                }else{
                                        r_degree1=ra1;
                                        r_degree=r_degree-ra1;
                                        r_degree2=min(r_degree, ra2);
                                }
                        }
                }

		if(min(l_degree1 , r_degree1) + min(l_degree2, r_degree2) >= min(l_degree1, r_degree2) + min( l_degree2, r_degree1)){
			if(min(l_degree1 , r_degree1)!=0){
			       new_AS_SV_i=nrow(new_AS_SV)+1;
				for(j in 1:15){
					new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
				}
				new_AS_SV[new_AS_SV_i,16]=min(l_degree1 , r_degree1);
                                new_AS_SV[new_AS_SV_i,18]=min(l_degree1 , r_degree1);
			}
                        if(min(l_degree2 , r_degree2)!=0){
                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i,17]=min(l_degree2 , r_degree2);
                                new_AS_SV[new_AS_SV_i,19]=min(l_degree2 , r_degree2);
                        }
		}else{
                       if(min(l_degree1 , r_degree2)!=0){
                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i,16]=min(l_degree1 , r_degree2);
                                new_AS_SV[new_AS_SV_i,19]=min(l_degree1 , r_degree2);
                        }
                        if(min(l_degree2 , r_degree1)!=0){
                               new_AS_SV_i=nrow(new_AS_SV)+1;
                                for(j in 1:15){
                                        new_AS_SV[new_AS_SV_i,j]=AS_SV[i,j];
                                }
                                new_AS_SV[new_AS_SV_i,17]=min(l_degree2 , r_degree1);
                                new_AS_SV[new_AS_SV_i,18]=min(l_degree2 , r_degree1)
                        }


		}		
	}
}

write.table(new_AS_SV,"SVs.AS_SV.haplotype_phased",quote=F, sep="\t", col.names=F, row.name=F)
