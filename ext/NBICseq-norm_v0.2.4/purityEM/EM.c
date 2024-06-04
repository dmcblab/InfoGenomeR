

#include "EM.h"

THETA new_THETA(int n_in, int EV, double an, double sx)
{	THETA theta = malloc(sizeof(struct THETA_t));
	int i, n;
	n = n_in;
	if(n<0) {fprintf(stderr,"Error in new_THETA, n (=%d) must be nonnegative\n",n);exit(1);}
	if(n<2) n = 2;
	if(theta==NULL){
		fprintf(stderr,"Error, memory allocation failed\n");
		exit(1);
		}
	theta->n = n;
	theta->p = (double *) malloc(sizeof(double)*(n+5)); /*make sure at least 5 elements*/
	theta->mu = (double *) malloc(sizeof(double)*(n+5));
	theta->mu_new = (double *) malloc(sizeof(double)*(n+5));
	theta->tauv = (double *) malloc(sizeof(double)*(n+5));
	theta->index_nonzero = (int *) malloc(sizeof(int)*(n+5));
	theta->psd = (double *) malloc(sizeof(double)*(n+5));
	theta->pvalue = (double *) malloc(sizeof(double)*(n+5));
	theta->size_index = 0;
	theta->len_index = 0;
	theta->WF = 0;

	if(theta->p==NULL || theta->mu == NULL || theta->mu_new == NULL || theta->tauv ==NULL || theta->index_nonzero==NULL || theta->psd==NULL || theta->pvalue==NULL){
		fprintf(stderr,"Error in new_THETA, memory allocation failed\n");
		exit(1);
		}

	theta->alpha = 0.0;
	theta->xi = 0.0;
	theta->tau = 0.01;
	theta->omega = 0.0;
	theta->D = 2;
	theta->ploidy = 2;
	theta->pi = 0;
	theta->num_param = theta->n+3;
	for(i=0;i<n+1;i++){
		theta->p[i] = 1.0/(n+1);
		theta->mu[i] = ((double) i)/2.0;
		theta->mu_new[i] = theta->mu[i] + theta->alpha*(1 - theta->mu[i]) + theta->xi;
		theta->tauv[i] = theta->tau;
		theta->index_nonzero[i] = -1;
		theta->psd[i] = 0.0001;
		theta->pvalue[i] = 0.0;
		}

	if(EV!=0) theta->EV = 1; else theta->EV = 0;

	theta->bic = 1e9;
	theta->lk = -1e9;

	theta->an = an;
	theta->sx = sx;
	if(theta->an<0.0) theta->an = 0.0;

	return theta;
}

void copy_THETA(THETA theta_in, THETA theta)
{	int i;

	if(theta_in->n != theta->n){
		fprintf(stderr,"Error in copy_THETA: theta1 and theta2 are of different dimension;\n");
		fprintf_THETA(theta_in, stderr,1);
		fprintf_THETA(theta, stderr,1);
		exit(1);
		}
	for(i=0;i<theta->n+1;i++){
		theta->p[i] = theta_in->p[i];
		theta->mu[i] = theta_in->mu[i];
		theta->mu_new[i] = theta_in->mu_new[i];
		theta->index_nonzero[i] = theta_in->index_nonzero[i];
		theta->tauv[i] = theta_in->tauv[i];
		theta->psd[i] = theta_in->psd[i];
		theta->pvalue[i] = theta_in->pvalue[i];
		}
	theta->alpha = theta_in->alpha;
	theta->xi = theta_in->xi;
	theta->tau = theta_in->tau;
	theta->omega = theta_in->omega;
	theta->D = theta_in->D;
	theta->ploidy = theta_in->ploidy;
	theta->pi = theta_in->pi;
	theta->len_index = theta_in->len_index;
	theta->size_index = theta_in->size_index;
	theta->num_param = theta_in->num_param;

	theta->bic = theta_in->bic;
	theta->lk = theta_in->lk;
	theta->EV = theta_in->EV;
	
	theta->an = theta_in->an;
	theta->sx = theta_in->sx;

	return;
}

THETA duplicate_THETA(THETA theta_in)
{	THETA theta=NULL;
	int i;

	if(theta_in==NULL) return theta;
	
	theta = new_THETA(theta_in->n, theta_in->EV, theta_in->an, theta_in->sx);

	for(i=0;i<theta->n+1;i++){
		theta->p[i] = theta_in->p[i];
		theta->mu[i] = theta_in->mu[i];
		theta->mu_new[i] = theta_in->mu_new[i];
		theta->index_nonzero[i] = theta_in->index_nonzero[i];
		theta->tauv[i] = theta_in->tauv[i];
		theta->psd[i] = theta_in->psd[i];
		theta->pvalue[i] = theta_in->pvalue[i];
		}
	theta->alpha = theta_in->alpha;
	theta->xi = theta_in->xi;
	theta->tau = theta_in->tau;
	theta->omega = theta_in->omega;
	theta->D = theta_in->D;
	theta->ploidy = theta_in->ploidy;
	theta->pi = theta_in->pi;
	theta->len_index = theta_in->len_index;
	theta->size_index = theta_in->size_index;
	theta->num_param = theta_in->num_param;

	theta->bic = theta_in->bic;
	theta->lk = theta_in->lk;
	theta->EV = theta_in->EV;
	
	theta->an = theta_in->an;
	theta->sx = theta_in->sx;
	return theta;
}

THETA del_component(THETA theta, int j)
{	THETA theta_new=NULL;
	int i;
	theta_new = new_THETA(theta->n, theta->EV, theta->an, theta->sx);

	
	copy_THETA(theta,theta_new);

	if(j>=0&&j<=theta->n && theta->n > 1){
		for(i=j;i<theta->n;i++){
			theta_new->p[i] = theta->p[i+1];
			theta_new->tauv[i] = theta->tauv[i+1];
			theta_new->psd[i] = theta->psd[i+1];
			theta_new->pvalue[i] = theta->pvalue[i+1];
			}
		theta_new->n--;
		}

	if(update_THETA(theta_new, theta_new->alpha, theta_new->xi, theta_new->tau, theta_new->p, theta_new->tauv)!=0){
		fprintf(stderr,"Deleting the %dth component\n",j);
		fprintf(stderr,"The original theta is\n");
		fprintf_THETA(theta,stderr,1);
		}
	compPloidy(theta_new);
	return theta_new;
}


THETA setzero_component(THETA theta, int j)
{       THETA theta_new=NULL;
        theta_new = new_THETA(theta->n, theta->EV, theta->an, theta->sx);


        copy_THETA(theta,theta_new);
	if(j<theta_new->n+1 && j>=0) theta_new->p[j] = 0.0;
        if(update_THETA(theta_new, theta_new->alpha, theta_new->xi, theta_new->tau, theta_new->p, theta_new->tauv)!=0){
		fprintf(stderr,"setting %dth component as zero\n",j);
		fprintf(stderr,"Original theta is\n");
		fprintf_THETA(theta,stderr,1);
		}
        compPloidy(theta_new);
        return theta_new;
}




void fprintf_THETA(THETA theta, FILE *out, int header)
{	int i;

	if(theta==NULL) return;
	compPloidy(theta);
	fprintf(out, "EV\tWF\tploidy\tpurity\tD\tomega\tbic\tloglk\talpha\txi\n");
	fprintf(out,"%d\t%d\t%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",theta->EV,theta->WF,theta->ploidy,1-theta->pi,theta->D,theta->omega,theta->bic,theta->lk,theta->alpha,theta->xi);
	for(i=0;i<theta->n+1;i++){
		fprintf(out,"tau%d\t",i);
		}
	fprintf(out,"\n");
	for(i=0;i<theta->n+1;i++){
		fprintf(out,"%g\t",theta->tauv[i]);
		}
	fprintf(out,"\n");

	for(i=0;i<theta->n+1;i++){
		fprintf(out,"p%d\t",i);
		}       
        fprintf(out,"\n");
	for(i=0;i<theta->n+1;i++){
		fprintf(out,"%g\t",theta->p[i]);
		}
	fprintf(out,"\n");

	if(theta->WF==1){
		for(i=0;i<theta->n+1;i++){
			fprintf(out,"%g\t",theta->psd[i]);
			}
		fprintf(out,"\n");
		for(i=0;i<theta->n+1;i++){
			fprintf(out,"%g\t",theta->pvalue[i]);
			}
		fprintf(out,"\n");
		}

	return;
}

int update_THETA(THETA theta, double alpha, double xi, double tau, double *p, double *tauv)
{	int i;
	double ss;
	const double EPS = 1e-10, tau_l = 1e-10;
	int rslt = 0; /*returned value, nonzero value means error*/

	if(tauv!=NULL && tau<=0.0) {fprintf(stderr,"Error in update_THETA: nonpositive tau value %g\n",tau); rslt = 1;}
	if(alpha < 1.0 - EPS) theta->alpha = alpha; else theta->alpha = 1.0;
	theta->xi = xi;
	theta->tau = tau;
	theta->len_index = 0;

	ss = 0.0;
	for(i=0;i<theta->n+1;i++){
		if(p[i]<0.0) {fprintf(stderr,"Error in update_THETA: negative probability p[%d]=%g\n",i,p[i]); rslt = 1; /*exit(1)*/}
		ss += p[i];
		theta->mu_new[i] = theta->mu[i] + theta->alpha*(theta->mu[2]-theta->mu[i]) + theta->xi;
		}
	if(ss<=0.0) {
		fprintf(stderr,"Error in update_THETA: zero total probability\n");
		ss = 1.0;
		rslt = 1;
		}
	if(theta->EV!=1 && tauv!=NULL){
		for(i=0;i<theta->n+1;i++){
			theta->tauv[i] = tauv[i];
			if(tauv[i] <= 0.0 && p[i]>0.0) {
				//fprintf(stderr, "Error in update_THETA: negative vairance value tauv[%d]=%g\n",i,tauv[i]); 
				theta->tauv[i] = tau_l;
				//rslt = 1;
				}
			}
		}else{
		for(i=0;i<theta->n+1;i++){
			theta->tauv[i] = tau;
			}
		}

	if(theta->alpha < 1.0 - EPS){
		for(i=0;i<theta->n+1;i++){
			theta->p[i] = p[i]/ss;
			if(theta->p[i]>0.0){
				theta->index_nonzero[theta->len_index] = i;
				theta->len_index ++;
				}
			}
		}else{
		if(theta->n < 2){
			theta->n = 2;
			}
		for(i=0;i<theta->n+1;i++){
			theta->p[i] = 0.0;
			//theta->tauv[i] = 0.0;
			}
		theta->p[2] = 1.0;
		//theta->tauv[2] = 1.0;
		//theta->tau = 1.0;
		}


	return rslt;
}

void compPloidy(THETA theta)
{	int i, id_maxp, dff;
	double tmp;

	theta->D = 0.0; theta->omega = 0.0;	
	id_maxp = 0;
	for(i=0;i<theta->n+1;i++){
		theta->D += 2*theta->p[i]*theta->mu[i];
		if(theta->p[i] > theta->p[id_maxp]) id_maxp = i;
		}
	theta->ploidy = id_maxp;
	theta->omega = 1 - (1- theta->alpha)*theta->D/2; /*this formula will guranntee omega <=1.0, but not necessarily >=0*/
	theta->pi = theta->omega*theta->D/(2 + (theta->D - 2)*theta->omega);

	if(theta->EV==1){
		if(theta->alpha == 1.0) theta->num_param = 2.0; else theta->num_param = theta->len_index+2.0;// + exp(-theta->pi); /*the term exp(-theta->pi) is to penalize the purity (=1-pi) away from 1*/
		}else{
		if(theta->alpha == 1.0) theta->num_param = 2.0; else theta->num_param = 2*theta->len_index + 2.0;// + exp(-theta->pi);
		}
	dff = theta->ploidy -2;
	if(dff<0) dff = -dff;

	theta->num_param += 2*dff;
	tmp = theta->pi;
	tmp = tmp > 0.0 ? tmp : -tmp; /*tmp = |1-purity|; this is used to penalize the model with small purity*/
	theta->num_param += 10*tmp;

	return;
}


void destroy_THETA(THETA theta)
{	if(theta==NULL) return;
	if(theta->mu!=NULL) {free(theta->mu); theta->mu=NULL;}
	if(theta->mu_new!=NULL) {free(theta->mu_new); theta->mu_new=NULL;}
	if(theta->p!=NULL) {free(theta->p); theta->p=NULL;}
	if(theta->tauv!=NULL) {free(theta->tauv); theta->tauv=NULL;}
	if(theta->index_nonzero!=NULL){free(theta->index_nonzero); theta->index_nonzero=NULL;}
	if(theta->psd!=NULL) {free(theta->psd); theta->psd=NULL;}
	if(theta->pvalue!=NULL){free(theta->pvalue); theta->pvalue=NULL;}

	free(theta); theta=NULL;
	return;
}


double dnorm(double x, double mu, double sigma2)
{	const double PI = 3.141592653589793238;
	double y;
	y = 1/sqrt(2*PI*sigma2)*exp(-1/(2*sigma2)*(x-mu)*(x-mu));
	return y;
}

double dnormMix(double x, double *mu, double *sigma2, double *p, int n)
{	double y;
	int i;
	y = 0.0;
	for(i=0;i<n;i++){
		if(p[i]>0.0 && sigma2[i]>0.0) y += dnorm(x, mu[i], sigma2[i])*p[i];
		}
	return y;
}



double loglk(double *x, int B, THETA theta)
{	double lk = 0.0, tmp;
	int b;
	for(b=0;b<B;b++){
		tmp = dnormMix(x[b],theta->mu_new,theta->tauv, theta->p, theta->n+1);
		lk += log(tmp);
		}

	return lk;
}

double regloglk(double *x, int B, THETA theta)
{	double rlk=0.0, pn, sx, psum;
	int i;

	sx = 0.0; psum = 0.0;
	for(i=0;i<theta->n+1;i++){
		//if(theta->p[i]>0.0 && theta->tauv[i]>0.0) {sx += theta->p[i]*theta->tauv[i]; psum += theta->p[i];}
		sx = sx > theta->tauv[i] ? sx : theta->tauv[i];
		}
	//if(psum>0.0) {sx = sx/psum;} else sx = theta->sx;

	rlk = loglk(x,B, theta);
	pn = 0.0;
	for(i=0;i<theta->n+1;i++){
		if(theta->p[i]>0.0 && theta->tauv[i]>0.0) pn += sx/theta->tauv[i] + log(theta->tauv[i]);
		}
	pn = -theta->an * pn;
	//fprintf(stderr,"theta->an = %g; theta->sx = %g; pn=%g\n",theta->an, theta->sx ,pn);
	rlk = rlk + pn;	
	return rlk;
	}

double fisher_pi(double *x, int B, THETA theta, int j)
{	int j0,b,i;
	double rslt, tmp1, tmp2;
	
	if(j<0 || j > theta->n) {fprintf(stderr,"Error in fisher_pi: j out of bound\n");exit(1);}
	if(theta->p[j]<=0 || theta->tauv[j]<=0) return 1;

	if(theta->p[theta->ploidy]>=1.0) return 1e10;
	if(j!=theta->ploidy) j0 = theta->ploidy; 
	else {
		j0 = 0;
		while(j0<theta->n+1) {
			if(j0!=theta->ploidy && theta->p[j0]>0.0 && theta->tauv[j0]>0.0) break;
			j0++;
			}
		if(j0==theta->n+1) return 1e10;
		}

	rslt = 0.0;
	for(b=0;b<B;b++){
		tmp1 = 0.0;
		for(i=0;i<theta->n+1;i++){
			if(theta->p[i]>0.0 && theta->tauv[i]>0.0) tmp1 += theta->p[i] * dnorm(x[b],theta->mu_new[i],theta->tauv[i]);
			}
		tmp1 = tmp1*tmp1;
		tmp2 = dnorm(x[b],theta->mu_new[j],theta->tauv[j]) - dnorm(x[b],theta->mu_new[j0],theta->tauv[j0]);
		tmp2 = tmp2*tmp2;
		rslt += tmp2/tmp1;
		}

	return rslt;

}

int assignTHETA_Pvalue(double *x, int B, THETA theta)
{	int i;
	double z;
	if(theta==NULL||B<0||x==NULL) return 1;
	theta->WF = 1;
	for(i=0;i<theta->n+1;i++){
		theta->psd[i] = 1.0/sqrt(fisher_pi(x,B,theta,i)+1e-6);
		z = theta->p[i]/(theta->psd[i]+1e-4); /*in this case the zscore must be >0*/
		if(z<100) theta->pvalue[i] = 2.0*pchisq(1,z*z,0); else theta->pvalue[i] = 1.6e-23;
		if(theta->pvalue[i]>1.0) theta->pvalue[i] = 1.0;
		}

	return 0;
}


void Erho(double *rho, double *x, int B, int n, THETA thetat)
{	int b, i,k;
	double rowSum;

	for(b=0;b<B; b++){
		rowSum = 0.0;
		//for(i=0;i<n+1;i++){
		for(k=0;k<thetat->len_index;k++){
			i = thetat->index_nonzero[k];
			rho[b*(n+1)+i] = thetat->p[i]*dnorm(x[b],thetat->mu_new[i],thetat->tauv[i]);
			rowSum += rho[b*(n+1)+i];
			}
		if(rowSum > 0.0){
			//for(i=0;i<n+1;i++){
			for(k=0;k<thetat->len_index;k++){
				i = thetat->index_nonzero[k];
				rho[b*(n+1)+i] = rho[b*(n+1)+i]/rowSum;
				}
			}
		}

	return;
}

int Mstep_EV(double *rho, double *x, int B, int n, THETA thetat, THETA thetanew, int alpha_fix)
{	double a_alpha, b_alpha, c_alpha, c_xi;
	double tmp, det;
	int b, i,k;
	const double EPS = 1e-10;
	
	a_alpha = 0.0;
	b_alpha = 0.0;
	c_alpha = 0.0;
	c_xi = 0.0;

	for(b=0;b<B;b++){
		for(k=0;k<thetat->len_index;k++){
			i = thetat->index_nonzero[k];
			tmp = rho[b*(n+1)+i]*(thetat->mu[2] - thetat->mu[i]);
			a_alpha += tmp*(thetat->mu[2] - thetat->mu[i]);
			b_alpha += tmp;
			c_alpha += tmp*(thetat->mu[i] - x[b]);
			c_xi += rho[b*(n+1)+i]*(thetat->mu[i] - x[b]);
			}
		}
	/*If thetat->alpha is effectively 1*/
	if(thetat->alpha - 1.0 > -EPS && thetat->alpha -1.0 < EPS){
		thetanew->xi = 0.0;
		thetanew->alpha = 1.0;
		thetanew->tau = 0.0;
		thetanew->EV = 1;
		for(b=0;b<B;b++){
			thetanew->xi += x[b];
			}
		thetanew->xi = thetanew->xi/B;
		for(b=0;b<B;b++){
			thetanew->tau += (x[b]-thetanew->xi)*(x[b]-thetanew->xi);
			}
		thetanew->tau = thetanew->tau/B;
		thetanew->xi = 1.0 - thetanew->xi;
		
		for(i=0;i<thetanew->n+1;i++){
			thetanew->p[i] = 0.0;
			}
		thetanew->p[2] = 1.0;
		if(update_THETA(thetanew, thetanew->alpha, thetanew->xi, thetanew->tau, thetanew->p, thetanew->tauv)!=0){
			fprintf(stderr,"Error in Mstep_EV: the starting theta value is\n");
			fprintf_THETA(thetat,stderr,1);
			//exit(1);
			}
		return 0; /*thetanew has been updated, exit this function*/
		}

	/*If thetat->alpha is not 1*/
	if(alpha_fix!=1){
		det = B*a_alpha - b_alpha*b_alpha;
		if(det==0.0) return 1; /*Error*/
		else {
			thetanew -> alpha = -(B*c_alpha-b_alpha*c_xi)/det;
			thetanew -> xi = - (a_alpha*c_xi - b_alpha*c_alpha)/det;
			}
		}else{
		thetanew -> alpha = thetat->alpha;
		thetanew -> xi = - (c_xi + b_alpha*thetat->alpha)/B;
		}
	
	//update_THETA(thetanew, thetanew->alpha, thetanew->xi, thetanew->tau, thetanew->p, thetanew->tauv);

	/*update thetanew->tau*/
	thetanew -> tau = 0.0;
	for(b=0; b<B; b++){
		//for(i=0;i<n+1;i++){
		for(k=0;k<thetat->len_index;k++){
			i = thetat->index_nonzero[k];
			thetanew->tau += rho[b*(n+1)+i]*(x[b] - thetanew->mu_new[i])*(x[b] - thetanew->mu_new[i]);
			}
		}
	thetanew -> tau = thetanew->tau/B;

	/*update thetanew->p*/
	for(i=0;i<n+1;i++){
		thetanew->p[i] = 0.0;
		}
	//for(i=0;i<n+1;i++){
	for(k=0;k<thetat->len_index;k++){
		i = thetat->index_nonzero[k];
		thetanew->p[i] = 0.0;
		for(b=0;b<B;b++){
			thetanew->p[i] += rho[b*(n+1)+i];
			}
		thetanew->p[i] = thetanew->p[i]/B;
		}

	if(update_THETA(thetanew, thetanew->alpha, thetanew->xi, thetanew->tau, thetanew->p, thetanew->tauv)!=0){
		fprintf(stderr,"Error in Mstep_EV: the starting theta value is\n");
		fprintf_THETA(thetat,stderr,1);
		//exit(1);
		}

	return 0;
}



int Mstep_UEV(double *rho, double *x, int B, int n, THETA thetat, THETA thetanew, int alpha_fix)
{	double a_alpha, b_alpha, c_alpha, b_xi, c_xi;
	double tmp, det;
	int b=0, i=0,k=0;
	static double *rho_sum = NULL;
	static int size_rho_sum = 0;

	if(size_rho_sum < thetat->n+10){
		size_rho_sum = thetat->n+20;
		rho_sum = (double *) realloc(rho_sum, sizeof(double)*size_rho_sum);
		if(rho_sum==NULL) {fprintf(stderr,"Error in Mstep: memory allocation failed\n");exit(1);}
		}

	copy_THETA(thetat,thetanew);

	a_alpha = 0.0;
	b_alpha = 0.0;
	c_alpha = 0.0;
	b_xi = 0.0;
	c_xi = 0.0;

	for(b=0;b<B;b++){
		for(k=0;k<thetat->len_index;k++){
			i = thetat->index_nonzero[k];
			tmp = rho[b*(n+1)+i]*(thetat->mu[2] - thetat->mu[i])/thetat->tauv[i];
			a_alpha += tmp*(thetat->mu[2] - thetat->mu[i]);
			b_alpha += tmp;
			c_alpha += tmp*(thetat->mu[i] - x[b]);
			c_xi += rho[b*(n+1)+i]*(thetat->mu[i] - x[b])/thetat->tauv[i];
			b_xi += rho[b*(n+1)+i]/thetat->tauv[i];
			}
		}
	/*this function will not consider the case thetat->alpha = 1.0*/

	/*If thetat->alpha is not 1*/
	if(alpha_fix!=1){
		det = b_xi*a_alpha - b_alpha*b_alpha;
		if(det==0.0) return 1; /*Error*/
		else {
			thetanew -> alpha = -(b_xi*c_alpha-b_alpha*c_xi)/det;
			thetanew -> xi = - (a_alpha*c_xi - b_alpha*c_alpha)/det;
			}
		}else{
		thetanew -> alpha = thetat->alpha;
		thetanew -> xi = - (c_xi + b_alpha*thetat->alpha)/b_xi;
		}

	if(update_THETA(thetanew, thetanew->alpha, thetanew->xi, thetanew->tau, thetanew->p, thetanew->tauv)!=0){
		fprintf(stderr,"Error in Mstep_UEV: the starting theta value is\n");
		fprintf_THETA(thetat,stderr,1);
		//exit(1);
		}

	/*update thetanew->tauv*/
	for(k=0;k<thetanew->n+1;k++){
		thetanew->tauv[k] = 0.0;
		rho_sum[k] = 0.0;
		}
	for(b=0;b<B;b++){
		for(k=0;k<thetat->len_index;k++){
			i = thetat->index_nonzero[k];
			thetanew->tauv[i] += rho[b*(n+1)+i]*(x[b] - thetanew->mu_new[i])*(x[b] - thetanew->mu_new[i]);
			rho_sum[i] += rho[b*(n+1)+i];
			}
		thetanew->tauv[i] += 2*thetanew->an*thetanew->sx;
		rho_sum[i] += 2*thetanew->an;
		}
	for(k=0;k<thetat->len_index;k++){
		i = thetat->index_nonzero[k];
		thetanew->tauv[i] = thetanew->tauv[i]/rho_sum[i];
		}

	/*update thetanew->p*/
	for(i=0;i<n+1;i++){
		thetanew->p[i] = 0.0;
		}
	//for(i=0;i<n+1;i++){
	for(k=0;k<thetat->len_index;k++){
		i = thetat->index_nonzero[k];
		thetanew->p[i] = 0.0;
		for(b=0;b<B;b++){
			thetanew->p[i] += rho[b*(n+1)+i];
			}
		thetanew->p[i] = thetanew->p[i]/B;
		}

	//update_THETA(thetanew, thetanew->alpha, thetanew->xi, thetanew->tau, thetanew->p);
	if(update_THETA(thetanew, thetanew->alpha, thetanew->xi, thetanew->tau, thetanew->p, thetanew->tauv)!=0){
		fprintf(stderr,"Error in Mstep_UEV: the starting theta value is\n");
		fprintf_THETA(thetat,stderr,1);
		//exit(1);
		}

	return 0;		
}

int Mstep(double *rho, double *x, int B, int n, THETA thetat, THETA thetanew, int alpha_fix)
{	const double EPS = 1e-10;
	int rslt;
	int iteration = 0;

	if(thetat->alpha - 1.0 > -EPS && thetat->alpha -1.0 < EPS){
		thetat->EV = 1;
		thetat->alpha = 1.0;
		}
	if(thetat->EV==1){
		rslt = Mstep_EV(rho,x,B,n, thetat,thetanew,alpha_fix);
		}else{
		THETA theta1=NULL, theta2=NULL;
		int flag = 0, i;
		double dff;
		rslt = Mstep_UEV(rho,x,B,n,thetat,thetanew,alpha_fix);
		if(iteration==1){
			theta1 = new_THETA(thetanew->n, thetanew->EV, thetanew->an, thetanew->sx);
			theta2 = new_THETA(thetanew->n, thetanew->EV, thetanew->an, thetanew->sx);
			copy_THETA(thetanew,theta1);
			copy_THETA(thetanew,theta2);
			i = 0;
			do{
				rslt = Mstep_UEV(rho,x,B,n,theta1,theta2,alpha_fix);
				dff = diff_theta(theta1,theta2);
				if(dff < 1e-6){flag = 1;}
				i++;
				copy_THETA(theta2,theta1);
				}while(flag==0 && i < 100);
			copy_THETA(theta1,thetanew);
			destroy_THETA(theta1); theta1= NULL;
			destroy_THETA(theta2); theta2= NULL;
			}
		}	

	return rslt;
}


double diff_theta(THETA theta1, THETA theta2)
{	double dff=0.0;
	int i;
	if(theta1->n != theta2->n || theta1->EV!=theta2->EV) {
		fprintf(stderr,"Error in diff_theta: theta1 and theta2 are not comparable\n");
		fprintf_THETA(theta1,stderr,1);
		fprintf_THETA(theta2,stderr,1);
		exit(1);
		}

	for(i=0;i<theta1->n+1;i++){
		dff = fmax(dff,fabs(theta1->p[i]-theta2->p[i]));
		}
	dff = fmax(dff, fabs(theta1->alpha - theta2->alpha));
	dff = fmax(dff, fabs(theta1->xi - theta2->xi));
	if(theta1->EV==1 && theta2->EV==1) {
		dff = fmax(dff, fabs(theta1->tau - theta2->tau));
		}else {  
		for(i=0;i<theta1->n+1;i++){
			dff = fmax(dff,fabs(theta1->tauv[i]-theta2->tauv[i]));
			}
		}

	return dff;
}


int EM_init(double *x, int B, THETA theta_ini, THETA theta_em, int max_it, double epsilon, int alpha_fix_in)
{	int k, flag, rslt;
	double dff;
	static double *rho=NULL;
	static int size = 0;
	THETA theta1, theta2;
	int alpha_fix;
	const double EPS = 1e-10;

	alpha_fix = alpha_fix_in;
	if(theta_ini->alpha>1.0){
		fprintf(stderr,"Error in EM_init: the initial value of alpha must be less than or equal to 1\n");
		return 1; /*error, return nonzero*/
		}
	theta1 = new_THETA(theta_ini->n, theta_ini->EV, theta_ini->an, theta_ini->sx);
	theta2 = new_THETA(theta_ini->n, theta_ini->EV, theta_ini->an, theta_ini->sx);
	copy_THETA(theta_ini, theta1);
	if(rho==NULL){
		size = B*(theta_ini->n+1)+10;
		rho = (double *)malloc(sizeof(double)*size);
		}else{
		if(size < B*(theta_ini->n+1)){
			size = B*(theta_ini->n+1)+10;
			rho = (double *) realloc(rho, sizeof(double)*size);
			}
		}
	if(rho==NULL){
		fprintf(stderr,"Error in EM_init: memory allocation failed\n");
		exit(1);
		}

	flag=0;
	k = 0;
	while(flag==0&&k<max_it){
		//if(theta1->alpha < EPS){theta1->alpha = 0.0; alpha_fix = 1;}
		if(theta1->alpha - 1.0 > -EPS){	
				theta1->alpha = 1.0; alpha_fix = 1; theta1->EV = 1; theta1->tau = 0.1;
				//update_THETA(theta1, theta1->alpha, theta1->xi, theta1->tau, theta1->p);
				if(update_THETA(theta1, theta1->alpha, theta1->xi, theta1->tau, theta1->p, theta1->tauv)!=0){
					fprintf(stderr,"Error in EM_ini: k =%d\ntheta_ini is\n",k);
					fprintf_THETA(theta_ini,stderr,1);
					fprintf(stderr,"\ntheta1 is\n");
					fprintf_THETA(theta1,stderr,1);
					//exit(1);
					}
				}
		Erho(rho, x, B, theta1->n, theta1);
		rslt = Mstep(rho,x, B, theta1->n, theta1, theta2, alpha_fix);
		//rslt = Mstep_EV(rho,x, B, theta1->n, theta1, theta2, alpha_fix);
		if(rslt!=0) return rslt; /*error in Mstep*/
		dff = diff_theta(theta1, theta2);
		if(dff<epsilon){flag = 1;}
		copy_THETA(theta2, theta1); /*theta2 is the new theta, copy it to theta1 for next iteration*/
		k++;
		}

	copy_THETA(theta1, theta_em);
	destroy_THETA(theta1);
	destroy_THETA(theta2);
	if(update_THETA(theta_em,theta_em->alpha,theta_em->xi, theta_em->tau, theta_em->p,theta_em->tauv)!=0){
		fprintf(stderr,"Error in EM_ini: theta_em is\n");
		fprintf_THETA(theta_em,stderr,1);
		}
	compPloidy(theta_em);
	theta_em->lk = regloglk(x, B,theta_em);  theta_em->bic = -2*theta_em->lk + log(B)*theta_em->num_param;

	return 0;
}


static int which_p_min(THETA theta, int left) /*left indicate whether it is on the left of the peak or not; return theta->ploidy if no element before or after the peak*/
{	int i, start, end, ind_min;

	compPloidy(theta);
	if(left==1){
		start = 0; end = theta->ploidy;
		}else{
		start = theta->ploidy;
		end = theta->n;
		}

	ind_min = start;
	for(i=start;i<=end;i++){
		if(theta->p[i] < theta->p[ind_min]) ind_min = i;
		}

	return ind_min;
}

/*static int which_p_min_nonzero(THETA theta, int left)
{	int i, ind_min, start, end;
	
	compPloidy(theta);
        if(left==1){
                start = 0; end = theta->ploidy;
                }else{
                start = theta->ploidy;
                end = theta->n;
                }


	ind_min = theta->ploidy;
	for(i=start;i<end;i++){
		if(theta->p[i] < theta->p[ind_min] && theta->p[i]>0.0) ind_min = i; 
		}

	return ind_min;
}*/

double peakHeight(THETA theta, int i)
{	double peak;
	if(theta==NULL || i>theta->n || i<0 || theta->p[i]<=0.0 || theta->tauv[i]<=0.0) peak = 0.0;
	else peak = theta->p[i]*dnorm(theta->mu_new[i],theta->mu_new[i],theta->tauv[i]);

	return peak;
}

int which_min_Peak(THETA theta,int left)
{	int i, start, end,i_min;
	double peak_min=0.0,peak=0.0;
	if(left==1) {start=0; end = theta->ploidy;}
	else {start=theta->ploidy; end = theta->n;}

	i_min = start;
	peak_min = peakHeight(theta,start);
	for(i=start+1;i<=end;i++){
		peak = peakHeight(theta,i);
		if(peak<peak_min) {peak_min = peak; i_min=i;}
		}
	return i_min;
}

THETA forceReduce(THETA theta_in, double *x, int B, int max_it, double epsilon)
{	int ind_min,alpha_fix=0;
	THETA theta1=NULL, theta2 = NULL;	

	ind_min = which_min_Peak(theta_in,1);
	if(theta_in->n>2 && ind_min!=theta_in->ploidy && theta_in->ploidy >2){
		//fprintf(stderr,"ind_min=%d\n",ind_min);
		//fprintf_THETA(theta_in,stderr,1);
		theta1 = del_component(theta_in, ind_min);
		//fprintf(stderr,"Deletion finished\n");
		compPloidy(theta1);
		theta2 = new_THETA(theta1->n, theta1->EV, theta1->an, theta1->sx);
		EM_init(x,B,theta1,theta2, max_it, epsilon, alpha_fix);
		//fprintf(stderr,"EM_init finished\n");
		destroy_THETA(theta1); theta1 = NULL;
		}else{
		theta2 = duplicate_THETA(theta_in);
		}


	return theta2;
}


void reduceModel(THETA theta_in, double *x, int B, int max_it, double epsilon)
{	//const double min_prob = 0.01;
	THETA theta1=NULL, theta2=NULL;// theta_ini;
	int alpha_fix = 0,flag=0, ind_min, left;
	double lk1, lk2, bic1, bic2;


	compPloidy(theta_in);
	ind_min=0;
	lk1 = regloglk(x, B,theta_in);
	bic1 = -2*lk1 + theta_in->num_param*log(B);
	flag = 0;
	for(left=0;left<=1;left++){
		flag = 0;
		ind_min = which_p_min(theta_in,left);
		if(left==1 && theta_in->ploidy<=2) flag = 1;
		if(theta_in->n<=2||ind_min==theta_in->ploidy) flag = 1;
		while(flag!=1){
			theta1 = del_component(theta_in, ind_min);
			compPloidy(theta1);
			theta2 = new_THETA(theta1->n, theta1->EV, theta1->an, theta1->sx);
			EM_init(x,B,theta1,theta2, max_it, epsilon, alpha_fix);
			lk2 = theta2->lk; bic2 = theta2->bic;
			if(bic2<bic1){
				theta_in->n--;
				copy_THETA(theta2,theta_in);
				bic1 = bic2;
				}else{
				flag = 1;
				}
			ind_min = which_p_min(theta_in,left);
			if(left==1 && theta_in->ploidy<=2) flag = 1;
			if(theta_in->n <=2 || ind_min==theta_in->ploidy) flag = 1; //theta->n should be at least 2
			destroy_THETA(theta1); theta1=NULL;
			destroy_THETA(theta2); theta2=NULL;
			}
		}

	/*
	flag = 0;
	while(flag!=1){
		ind_min = which_p_min_nonzero(theta_in,0); //only search the right side
		if(ind_min!=theta_in->ploidy){
			theta1 = setzero_component(theta_in, ind_min);
			compPloidy(theta1);
			theta2 = new_THETA(theta1->n, theta1->EV, theta1->an, theta1->sx);
			EM_init(x,B,theta1,theta2, max_it, epsilon, alpha_fix);
			if(theta2->bic < theta_in->bic){
				copy_THETA(theta2,theta_in);
				}else{
				flag = 1;
				}
			destroy_THETA(theta1); theta1=NULL;
			destroy_THETA(theta2); theta2=NULL;
			}else{
			flag = 1;
			}
		}*/

	return;
}





THETA EM_purity(double *x, int B, int max_it, double epsilon)
{	int n, i, b, k, flag;
	double *cnt, mu, tmp;
	THETA theta_ini, theta_max, theta_em, theta_bic; /*theta_bic will store the model chosen by bic criterion*/
	double alpha0 = 0.1, xi0 = 0.0, tau0=0.01, alpha_fix=0, sd0=0.1;
	double lk=0.0, lk_max=-1e10, bic=1e10, bic_min=1e10,lk_bic=-1e10;
	const double EPS = 1e-10;

	n = (int) 2.0*x[B-1] + 1;
	theta_ini = new_THETA(n,1,0,1);
	theta_max = new_THETA(n,1,0,1);
	theta_em = new_THETA(n,1,0,1);
	theta_bic = new_THETA(n,1,0,1);

	cnt = (double *) malloc(sizeof(double)*(n+2));
	for(i=0;i<n+2;i++){
		cnt[i]=0.0;
		}
	k=0; mu = 0.0;
	b = 0;
	for(k=0;k<n+2;k++){
		mu = ((double) k + 1.0)/2;
		while(b<B && x[b] < mu){
			cnt[k] += 1.0;
			b++;
			}
		}
	tmp = 0.0;
	for(i=0;i<theta_ini->n+1;i++){
		theta_ini->p[i] = cnt[i]+cnt[i+1];
		tmp += theta_ini->p[i];
		}
	if(tmp<=0.0){fprintf(stderr,"Error in EM: zero totol proability (%g) generated\n", tmp);exit(1);}
	for(i=0;i<theta_ini->n+1;i++){
		theta_ini->p[i] = theta_ini->p[i]/tmp;
		}
	flag = 0;
	for(alpha0=0.0; alpha0 < 1.0 + EPS; alpha0 += 0.10){
		int num_param;
		if(alpha0 - 1.0 > -EPS || alpha0 < EPS) alpha_fix = 1; else alpha_fix = 0;
		if(alpha0 - 1.0 > -EPS) {num_param = 2; alpha0 = 1.0;}
		else if(alpha0 < EPS) {num_param = theta_ini->n + 2; alpha0 = 0.0;}
		else num_param = theta_ini->n + 3;

		for(sd0=0.1; sd0<=1.0; sd0+=0.1){
			tau0 = sd0*sd0;
			update_THETA(theta_ini, alpha0, xi0, tau0, theta_ini->p,theta_ini->tauv);
			EM_init(x, B,theta_ini,theta_em,max_it, epsilon, alpha_fix);
			lk = regloglk(x, B,theta_em); bic = -2*lk + log(B)*num_param;
			if(flag==0){lk_max = lk; copy_THETA(theta_em, theta_max); flag=1; bic_min = bic; lk_bic = lk;}
			else{
				if(lk>lk_max){lk_max = lk; copy_THETA(theta_em, theta_max);}
				if(bic < bic_min){bic_min = bic; lk_bic = lk; copy_THETA(theta_em, theta_bic);}
				}
			}
		}
	compPloidy(theta_bic);
	theta_bic->lk = lk_bic; theta_bic->bic = bic_min;
	//fprintf(stdout,"log likelihood: %g\n",lk_max);
	//fprintf_THETA(theta_max, stdout, 0);
	//fprintf(stdout,"minium bic: %g; Corresponding log likelihood is %g\n",bic_min,lk_bic);
	//fprintf_THETA(theta_bic, stdout, 0);
	//destroy_THETA(theta_ini);
	//destroy_THETA(theta_em);
	//destroy_THETA(theta_max);
	
	return theta_bic;
}




THETA_VEC new_THETA_VEC(int n)
{	THETA_VEC theta_v = malloc(sizeof(struct THETA_VEC_t));
	int i;

	if(n<1) theta_v->size = 1; else theta_v->size = n;

	theta_v->theta = (THETA *) malloc(sizeof(THETA)*theta_v->size);
	if(theta_v->theta==NULL){
		fprintf(stderr,"Error in new_THETA_VEC: memoery allocation failed\n");
		exit(1);
		}

	for(i=0;i<theta_v->size;i++){
		theta_v->theta[i] = NULL;
		}
	theta_v->lst = -1;

	return theta_v;
}

void destroy_THETA_VEC(THETA_VEC theta_v)
{	int i;
	if(theta_v==NULL) return;

	for(i=0;i<theta_v->size;i++){
		destroy_THETA(theta_v->theta[i]);
		theta_v->theta[i] = NULL;
		}
	theta_v->lst = 0;
	theta_v->size = 0;
	free(theta_v);
	theta_v = NULL;
	return;	
}

void extend_THETA_VEC(THETA_VEC theta_v, int n)
{	int i;

	if(theta_v==NULL){
		fprintf(stderr,"Error in enxtend_THETA_VEC: theta_v is NULL\n");
		exit(1);
		}

	if(n>theta_v->size){
		theta_v->theta = (THETA *) realloc(theta_v->theta, sizeof(THETA)*n);
		if(theta_v->theta==NULL){
			fprintf(stderr,"Error in extend_THETA_VEC: memoery allocation failed\n");
			exit(1);
			}
		for(i=theta_v->size;i<n;i++){
			theta_v->theta[i] = NULL;
			}
		theta_v->size = n;
		}

	return;
}


void set_value_THETA_VEC(THETA_VEC theta_v, THETA theta, int i)
{
	if(i<0){
		fprintf(stderr,"Error in set_value_THETA_VEC: i=%d is negative\n",i);
		exit(1);
		}

	if(theta_v==NULL){
		fprintf(stderr,"Error in set_value_THETA_VEC: theta_v is NULL\n");
		exit(1);
		}

	if(i>=theta_v->size){
		extend_THETA_VEC(theta_v, i+10);
		}
	if(theta_v->theta[i]!=NULL){
		destroy_THETA(theta_v->theta[i]); 
		theta_v->theta[i]=NULL;
		}
	theta_v->theta[i] = duplicate_THETA(theta);
	
	if(i>theta_v->lst && theta!=NULL) theta_v->lst = i;
	if(i==theta_v->lst && theta==NULL) theta_v->lst--;

	return; 	
}


int predict_copynumber(THETA theta, double x, double min_density)
{	static double *density = NULL;
	static int size_density = 0;
	int i;
	int index_max1, index_max2, sum;
	double pst;

	if(density==NULL||size_density<theta->n){
		size_density = 10 > theta->n+5 ? 10: theta->n+2;
		density = (double *) realloc(density,size_density*sizeof(double));
		if(density==NULL){
			fprintf(stderr,"Error in predict_copynumber: memory allocation failed\n");
			exit(1);
			}
		}

	index_max1 = 0;
	index_max2 = 0;
	sum = 0.0;
	for(i=0;i<theta->n+1;i++){
		if(theta->tauv[i] <=0.0){
			density[i] = 0.0;
			}else{
			density[i] = dnorm(x, theta->mu_new[i],theta->tauv[i]);
			}
		if(density[i] > density[index_max1]){
			index_max1 = i;
			}
		if(theta->p[i]*density[i] > theta->p[index_max2]*density[index_max2]){
			index_max2 = i;
			}
		sum += theta->p[i]*density[i];
		}

	if(density[index_max1] < min_density || density[index_max2] < min_density){ /*density too small, return -1*/
		return -1; 
		}

	pst = theta->p[index_max2]*density[index_max2]/sum;
	if(pst < min_density){
		return -1;
		}

	return index_max2;
}
