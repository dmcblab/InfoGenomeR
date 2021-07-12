args=commandArgs(T);

search_length=as.integer(args[4]);

#library(plot3D)
t=read.table(args[1]);
t$V2=abs(t$V2);
t$V4=abs(t$V4);
t=t[t$V10>0.8,]
if(nrow(t)==0){
	print(0);
}else{

	x=t[,2]-search_length;
	x=x%/%2000+search_length/2000;
	y=t[,4]-search_length;
	y=y%/%2000+search_length/2000;
	z=matrix(0,nrow=2*search_length/2000, ncol=2*search_length/2000)
	for(z_i in 1:length(x))
			z[x[z_i],y[z_i]]=z[x[z_i],y[z_i]]+t[z_i,10];

	z=z[z!=0];
	options(scipen = 999);
	if(length(z)>0 && max(round(z))>2){
		if(max(round(z))>4){
                        print(max(ppois(round(z),as.integer(args[2])*as.integer(args[3]))));
#			print(1);
		}else{
			print(max(ppois(round(z),as.integer(args[2])*as.integer(args[3]))));
		}
	}else{
		print(0);
	}
}
#png("test.png")
#hist3D(z=z, border="black", zlim=c(0,max(z)))
#dev.off()
