/*
Written by Ruibin Xi
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include "../lib/read.h"
#include <getopt.h>
#include <math.h>
#include <unistd.h>
#include "EM.h"
#include "gamma.h"

typedef struct args_purityEM_t{
        char *infile, *outfile;
        int max_it; /*maximum number of iterations*/
        double epsilon; /*when to stop the EM iterations*/
	double mploidy; /*a rough upper bound of the ploidy number*/
	int maxComp; /*number of components in the mixture model*/
	int qck; /*whether or not use maxComp method to estimate the parameters*/
	int nRS; /*number of random starting values for each possible value (for maxComp method)*/
	int subsample;
	int printV;
        FILE *in,*out;
        } args_purityEM; /*the argument passed to the program*/
static const int Nsubsample = 10000;

static args_purityEM option_assign(int argc, char **argv);
static void explain_command(char **argv);
static int cmpdouble(const void *a, const void *b);

static double MAD(double *x, int B);/*calculate the MAD of x, where x should have been sorted*/
static void removeOutlier(double *copyratio, int B, int *k1, int *k2); 
/*given B copy ratios (sorted nondecreasingly), 
  calculate k1 and k2 (k1 < k2) such that the first k1 elements and the last B-k2 elelments are likely outliers
*/

static double *filterOutlier(double *copyratio, int n, int *n_remain, THETA theta, double min_density, double min_comp_size);
/* given a model theta,
   filter the likely outliers, return a vector of remaining copy ratios.
   n_remain is the length of the remaining vector
   min_comp_size: the minimum component probability
*/

static THETA selectBestTheta(double *x, int B, THETA_VEC vct, THETA theta_bic); /*select the best model*/


static THETA create_init(double *x, int B, double omega, double D, double tau, int EV, double an, double sx); 
/*omega and D are defined as in the manuscript (D is the ploidy number); 
  x should be ordered nondecreasingly
   B is the length  of x
  EV: whether or not the variances should be the same
*/

static THETA create_init_new(double *x, int B, int numComp, int EV, double an, double sx);
/* x should be ordered nondecreasingly
   B is the length of x
   EV: whether or not the variances should be the same
   numComp: number of components in the mixture model
*/
static void update_bic_solution(THETA *ptheta_bic, THETA thetanew);
static void update_bic_solution_vec(THETA_VEC vec, THETA thetanew);

static THETA randomGen_init(int numComp, int EV, double an, double sx); /*randomly generate a starting value*/

static THETA EM_auto(double *x, int B, int max_it, double epsilon, double max_ploidy, double an, double sx, THETA_VEC *vct);

static THETA EM_auto_new(double *x, int B, int max_it, double epsilon, int maxComp, double an, double sx, THETA_VEC *vct);
/*the ith element of vct will be the best model corresponds to ploidy number i*/

static THETA EM_auto_randomstart(double *x, int B, int max_it, double epsilon, int maxComp, double an, double sx, int nRS, THETA_VEC *vct);
/*nRS: numver of random starting values*/

int main(int argc, char **argv)
{	args_purityEM args;
	double *copyratio, *x, *x1;
	int nrow, ncol, B, start, end, B1;
	double an, sx;
	THETA_VEC vct = NULL, vct1 = NULL;
	THETA theta_bic = NULL, theta1 = NULL, theta2 = NULL, theta_b=NULL;
	//int i;

	

	args = option_assign(argc,argv);
	copyratio = read_table(args.in, &nrow, &ncol, -1, 0);
	if(ncol!=1){
		fprintf(stderr,"<input> must be a one column file\n");
		exit(1);
		}
	if(nrow<10){
		fprintf(stderr,"Too few observations to make meaningfule inference\n");
		exit(2);
		}
	qsort(copyratio,nrow,sizeof(double),cmpdouble);
	if(copyratio[0]<0){
		fprintf(stderr,"the copy ratios must be nonnegative\n");
		exit(1);
		}
	//need a few lines to handle the special case of all of copyratio are zero
	removeOutlier(copyratio,nrow, &start, &end);
	B = end - start + 1;
	x = &copyratio[start];

	if(args.subsample==1 && B > Nsubsample){
		int seed;
		seed = genSeed();
		seed_set(seed);
		db_shuffle(x,B);
		B = Nsubsample;
		qsort(x,B,sizeof(double),cmpdouble);
		//fprintf(stderr,"B=%d\n",B);
		}
	
	vct = new_THETA_VEC(10);
	an = 1.0/sqrt(B+0.0);
	sx = MAD(x,B);
	if(args.qck == 1){
		theta1 = EM_auto_new(x, B, args.max_it, args.epsilon,args.maxComp,an,sx, &vct);
		update_bic_solution(&theta_bic, theta1);
		theta2 = EM_auto_randomstart(x, B, args.max_it, args.epsilon,args.maxComp,an,sx, args.nRS, &vct);
		update_bic_solution(&theta_bic, theta2);
		}else{
		theta_bic = EM_auto(x, B, args.max_it, args.epsilon,args.mploidy,an,sx, &vct);
		}

	/*
	if(args.printV==1){
		for(i=0;i<=vct->lst;i++){
			assignTHETA_Pvalue(x,B,vct->theta[i]);
			fprintf_THETA(vct->theta[i], args.out, 1);
			fprintf(args.out,"\n");
			}
		}*/

	theta_b = selectBestTheta(x,B, vct, theta_bic);
	//fprintf(stderr,"The best selected model is\n");
	fprintf_THETA(theta_b,args.out,1);

	/*filter the outliers and run the algorithm again*/
	x1 = filterOutlier(x, B, &B1, theta_b, 1e-6, 0.005);
	//for(i=0;i<B1;i++){
	//	fprintf(stdout,"%g\n",x1[i]);
	//	}

	
	
	if(theta1!=NULL) {destroy_THETA(theta1); theta1=NULL;}
	if(theta2!=NULL) {destroy_THETA(theta2); theta2=NULL;}
	if(theta_bic!=NULL) {destroy_THETA(theta_bic); theta_bic = NULL;}
	if(theta_b!=NULL) {destroy_THETA(theta_b); theta_b = NULL;}


	vct1 = new_THETA_VEC(10);
	theta_bic = NULL;
	theta1 = EM_auto_new(x1, B1, args.max_it, args.epsilon,args.maxComp,an,sx, &vct1);
	update_bic_solution(&theta_bic, theta1);
	theta2 = EM_auto_randomstart(x1, B1, args.max_it, args.epsilon,args.maxComp,an,sx, args.nRS, &vct1);
	update_bic_solution(&theta_bic, theta2);

	theta_b = selectBestTheta(x1,B1, vct1, theta_bic);
	fprintf_THETA(theta_b,args.out,1);

	return 0;
}


static args_purityEM option_assign(int argc, char **argv)
{	args_purityEM args;
	int c,option_index;
	int flag_mploidy=0, flag_maxComp=0;
	static struct option long_options[] =
		{ {"max_it", required_argument, 0, 0},
		  {"epsilon",required_argument, 0, 1},
		  {"mploidy",required_argument, 0, 2},
		  {"maxComp",required_argument, 0, 3},
		  {"subsample",no_argument,0,4},
		  {"nRS", required_argument, 0, 5},
		  {0,0,0,0}
		};
	args.max_it = 100;
	args.epsilon = 1e-6;
	args.mploidy = 4;
	args.maxComp = 8;
	args.qck = 1;
	args.nRS = 200;
	args.out = stdout;
	args.in = NULL;
	args.outfile = NULL;
	args.infile = NULL;
	args.subsample = 0;
	args.printV = 0;

	while((c=getopt_long(argc,argv,"po:h",long_options,&option_index))!=-1){
		switch(c){
			case 0:
				args.max_it = atoi(optarg);
				if(args.max_it < 10) args.max_it = 10;
				break;
			case 1:
				args.epsilon = atof(optarg);
				if(args.epsilon <= 0.0) args.epsilon = 1e-6;
				break;
			case 2:
				args.mploidy = atof(optarg);
				if(args.mploidy<2.0) args.mploidy = 2.0;
				flag_mploidy = 1;
				break;
			case 3:
				args.maxComp = atoi(optarg);
				if(args.maxComp<3) args.maxComp = 3;
				flag_maxComp = 1;
				break;
			case 4:
				args.subsample = 1;
				break;
			case 5: 
				args.nRS = atoi(optarg);
				if(args.nRS<1) args.nRS = 0;
				break;
			case 'o':
				args.outfile = strdup(optarg);
				break;
			case 'p':
				args.printV = 1;
				break;
			case 'h':
				explain_command(argv);
				exit(0);
				break;
			case '?': /* getopt_long already printed an error message. */
				exit(1);
				break;
			default: 
				abort();
			}
		}
	if(flag_mploidy==1 && flag_maxComp !=1){
		args.qck = 0;
		}


	if (argc - optind!=1){
		explain_command(argv);
		exit(1);
		}
	args.infile = strdup(argv[optind]);

	args.in = fopen(args.infile,"r");
	if(args.in==NULL){
		fprintf(stderr,"fopen: %s: %s\n", args.infile,strerror(errno));
		exit(1);
		}
	if(args.outfile!=NULL){
		args.out = fopen(args.outfile,"w");
		if(args.out==NULL){
			fprintf(stderr,"fopen: %s: %s\n", args.outfile,strerror(errno));
			exit(1);
			}
		}
	return args;
}


static void explain_command(char **argv)
{	fprintf(stderr,"Usage: %s [options] <input>\n",argv[0]);
	fprintf(stderr,"Options:\n");
	fprintf(stderr,"        -h : print this message\n");
	fprintf(stderr,"        -o <string>: the output name; default <stdout>\n");
	fprintf(stderr,"        -p : print the best model for each ploidy number\n");
	fprintf(stderr,"        --max_it <int>: maximum number of EM iterations; default 100\n");
	fprintf(stderr,"        --epsilon <float>: the convergence criterion; default 1e-6\n");
	fprintf(stderr,"        --mploidy: specify a rough upper bound of the ploidy number.\n"); 
	fprintf(stderr,"                   A smaller value (e.g. 4) is suggested as a larger value will make the program quite slow\n");
	fprintf(stderr,"        --maxComp: specify the maximum number of component in the mixture of normal model; Default 15\n");
	fprintf(stderr,"                   this option will overload the option --mploidy; this is the default method\n");
	fprintf(stderr,"        --nRS: the number of random starting values for the EM algorithm\n");
	fprintf(stderr,"        --subsample: if specified, subsample %d observations (if there are more than %d observations) for parameter estimation.\n", Nsubsample, Nsubsample);
	return;
}


static int cmpdouble(const void *a, const void *b)
        { double tmp1,tmp2;
          tmp1 = *((const double *) a);
          tmp2 = *((const double *) b);
          if(tmp1<tmp2) return -1;
          else if(tmp1>tmp2) return 1;
          else return 0;
        }

static double MAD(double *x, int B)
{	double *dff=NULL, mdn, mad;
	int size,i;
	const double scale=1.4826;

	if(B<0) {fprintf(stderr,"Error in MAD: B = %d is nongative\n",B);exit(1);}
	size = B+10;
	dff = (double *) malloc(sizeof(x)*(size));
	if(dff==NULL) {fprintf(stderr,"Error in MAD: memory allocation failed\n"); exit(1);}
	mdn = x[B/2];
	for(i=0;i<B;i++){
		dff[i] = x[i] - mdn;
		dff[i] = dff[i] >=0.0 ? dff[i] : -dff[i];
		}
	qsort(dff,B,sizeof(double),cmpdouble);

	mad = dff[B/2];
	if(dff!=NULL) {free(dff); dff = NULL;}

	return scale*mad;
}

static void removeOutlier(double *x, int B, int *k1, int *k2)
{	double median, mad, ub,lb;
	int i;
	median = x[B/2];
	mad = MAD(x,B);
	ub = median + 10*mad;
	lb = median - 10*mad;
	fprintf(stderr,"median = %g; MAD=%g, LB =%g, UB=%g, min=%g, max = %g\n",median,mad,lb, ub,x[0], x[B-1]);

	i = B-1;
	while(i>=B/2&&x[i] > ub){
		i--;
		}
	*k2 = i;
	i = 0;
	while(i<=B/2&&x[i]<lb){
		i++;
		}
	*k1 = i;

	return;
}

static THETA selectBestTheta(double *x, int B, THETA_VEC vec, THETA theta_bic)
{	int i, ploidy,flag1 = 0, flag2 = 0, start=0;
	const double minPeakPerc = 0.005, siglevel = 0.01;
	double peak, ploidyPeak, PeakThreshold;
	THETA theta_b=NULL;
	if(vec==NULL||vec->lst<=2) return theta_bic;

	start = 2;
	for(i=0;i<=vec->lst;i++){
		assignTHETA_Pvalue(x,B,vec->theta[i]);
		//if(vct->theta[start]==NULL && i>=2 && vct->theta[i]!=NULL && i<=theta_bic->ploidy) start = i;
		}
	//if(vct->theta[start]==NULL || start==theta_bic->ploidy) return theta_bic;
	assignTHETA_Pvalue(x,B,theta_bic);

	ploidy = theta_bic->ploidy;
	theta_b = duplicate_THETA(theta_bic);
	flag1 = 0; 
	while(flag1==0 && ploidy > 2){
		flag2 = 0; i = 0;
		ploidyPeak = theta_b->p[ploidy]*dnorm(theta_b->mu_new[ploidy],theta_b->mu_new[ploidy],theta_b->tauv[ploidy]);
		PeakThreshold = minPeakPerc*ploidyPeak;
		//fprintf(stderr,"The best model so far is\n");
		//fprintf_THETA(theta_b, stderr,1);
		//fprintf(stderr,"Ploidy peak is %g; minimum peak is %g\n",ploidyPeak, PeakThreshold);
		while(flag2==0&&i<theta_b->ploidy){ /*if any p of theta_b is too small or not significant*/
			peak = theta_b->p[i]*dnorm(theta_b->mu_new[i],theta_b->mu_new[i],theta_b->tauv[i]);
			//fprintf(stderr,"peak[%d]=%g\n",i,peak);
			if(peak <= PeakThreshold || theta_b->pvalue[i] >= siglevel) flag2 = 1;
			i++;
			}
		//fprintf(stderr,"flag2=%d\n",flag2);
		//fprintf(stderr,"\n");

		
		if(flag2 == 1){/*some of p is too small or not significant, look for the next best one from 2 to ploidy*/
			int ind_best = start;
			destroy_THETA(theta_b); theta_b=NULL;
			//theta_b = vec->theta[start];
			for(i=start;i<ploidy;i++){
				if(vec->theta[i]!=NULL && vec->theta[i]->bic <= vec->theta[ind_best]->bic) ind_best = i; //theta_b = vec->theta[i];
				}
			theta_b = duplicate_THETA(vec->theta[ind_best]);
			if(theta_b==NULL||theta_b->ploidy == ploidy) flag1 = 1; else ploidy = theta_b->ploidy;
			}else{
			flag1 = 1;
			}
		}

	return theta_b;

}


static THETA create_init(double *x, int B, double omega, double D, double tau, int EV, double an, double sx)
{	THETA theta_init=NULL;
	static int n, *cnt=NULL, size=0;
	int i,b;
	double dist;
	double EPS = 1e-10;
	if(D<=EPS){
		fprintf(stderr,"The ploidy number must be positive\n");
		exit(1);
		}
	if(tau<=EPS){
		fprintf(stderr,"Error in create_init: tau must be positive\n");
		exit(1);
		}
	if(omega < 1.0-EPS){
		dist = (1-omega)/D;
		n = (int) x[B-1]/dist+1.0;
		//if(n<2) n = 2;
		if(n>size){
			size = n+10;
			cnt = (int *) realloc(cnt,sizeof(int)*size);
			if(cnt==NULL){
				fprintf(stderr,"Error in create_init: memory allocation failed\n");
				exit(1);
				}
			}
			/*initialize cnt*/
			for(i=0;i<=n;i++){
				cnt[i] = 0;
				}
			/*cnt is the number of observations in [i*dist, (i+1)*dist)*/
			b = 0;
			for(i=0;i<=n;i++){
				while(b < B && x[b]<(i+1)*dist){
					cnt[i]++;
					b++;
					}
				}
			theta_init = new_THETA(n, EV, an, sx);
			theta_init->p[0] = (double) cnt[0];
			for(i=1;i<=n;i++){
				theta_init->p[i] = (double) (cnt[i-1]+cnt[i]);
				//fprintf(stderr,"%g\t",theta_init->p[i]);
				}
			//fprintf(stderr,"\n");
			theta_init->alpha = 1-2*(1-omega)/D;
			theta_init->xi = omega - theta_init->alpha;
			theta_init->tau = tau;
			//fprintf(stderr,"dist=%g, omega=%g, tau = %g, D=%g, n = %d\n",dist,omega, tau, D, n);
			while(theta_init->p[n] <=0.0 && theta_init->n > 2){ /*p[n] is zero, it will always be zero in the EM iterations, so no need to keep it*/
				theta_init->n --;
				}
			update_THETA(theta_init, theta_init->alpha, theta_init->xi, theta_init->tau, theta_init->p, theta_init->tauv);
		}else{
		n = 2;
		theta_init = new_THETA(n, EV, an, sx);
		theta_init->p[2] = 1.0; theta_init->p[0] = theta_init->p[1] = 0.0;
		theta_init->alpha = 1.0;
		theta_init->xi = 0.0;
		theta_init->tau = tau;
		update_THETA(theta_init, theta_init->alpha, theta_init->xi, theta_init->tau, theta_init->p, theta_init->tauv);
		}
	return theta_init;
}

static THETA create_init_new(double *x, int B, int numComp_in, int EV, double an, double sx)
{	THETA theta = NULL;
	int nc, i,b;
	static int *cnt=NULL, size=0;
	const double EPS = 1e-10;
	double dist;
	if(numComp_in <= 1) nc = 1; else nc = numComp_in;

	if(nc==1){
		theta = new_THETA(2, 1, an, sx);
		theta->alpha = 1.0;
		theta->p[0] = theta->p[1] = 0.0;
		theta->p[2] = 1.0;
		theta->tau = 0.01;
		theta->xi = 0.0;
		update_THETA(theta,theta->alpha,theta->xi,theta->tau,theta->p, NULL);
		return theta;
		}

	theta = new_THETA(nc-1, EV, an, sx);

	theta->alpha = 1-2*x[B-1]/(nc+1);
	if(theta->alpha <= 0.0) theta->alpha = 0.1; else if(theta->alpha > 1.0 - EPS) theta->alpha = 0.9;
	theta->xi = 0.0;
	dist = x[B-1]/nc;

 	if(nc+10>size){
		size = nc + 10;
		cnt = (int *) realloc(cnt, sizeof(int)*size);
		if(cnt==NULL){
			fprintf(stderr,"Error in create_init_new: memoery allocation failed\n");
			exit(1);
			}
		}
	/*initialize cnt*/
	for(i=0;i<=nc+1;i++){
		cnt[i] = 0;
		}
	/*cnt is the number of observations in [i*dist, (i+1)*dist)*/
	b = 0;
	for(i=0;i<=nc;i++){
		while(b < B && x[b]<(i+1)*dist){
			cnt[i]++;
			b++;
			}
		}
	theta->p[0] = (double) cnt[0];
	for(i=1;i<=nc-1;i++){
		theta->p[i] = (double) cnt[i-1] + cnt[i];
		}
	for(i=0;i<theta->n+1;i++){
		theta->tauv[i] = 0.1*0.1;
		}
	theta->tau = 0.1*0.1;
	update_THETA(theta,theta->alpha,theta->xi, theta->tau, theta->p, theta->tauv);

	return theta;
	
}

static THETA randomGen_init(int numComp, int EV, double an, double sx)
{	THETA theta = NULL;
	int nc,i;
	static int seed=5, flag=0;
	double *dalpha = NULL; /*parameter for Dirichelet distribution*/
	
	if(flag!=1){
		seed = genSeed();
		seed_set(seed);
		flag = 1;
		}

	if(numComp <= 1) nc = 1; else nc = numComp;
	if(nc==1){
                theta = new_THETA(2, 1, an, sx);
                theta->alpha = 1.0;
                theta->p[0] = theta->p[1] = 0.0;
                theta->p[2] = 1.0;
                theta->tau = 0.01;
                theta->xi = 0.0;
                update_THETA(theta,theta->alpha,theta->xi,theta->tau,theta->p, NULL);
                return theta;
                }

	theta = new_THETA(nc-1, EV, an, sx);

	theta->alpha = rand_lp();
	theta->xi = -0.5 + rand_lp();
	theta->tau = rgamma1(0.5,0.5);
	if(theta->tau< 0.0001) theta->tau = 0.0001;
	for(i=0;i<theta->n+1;i++){
		theta->tauv[i] = rgamma1(0.5,0.5);
		if(theta->tauv[i] < 0.0001) theta->tauv[i] = 0.0001;
		}

	dalpha = (double *) malloc(sizeof(double)*(theta->n+10));
	if(dalpha==NULL) {fprintf(stderr,"Error in randomGen_init: memory allocation failed\n"); exit(1);}
	for(i=0;i<theta->n+1;i++){
		dalpha[i] = 1.0;
		}
	rDirichlet(dalpha, theta->p, theta->n+1);
	
	if(update_THETA(theta,theta->alpha,theta->xi, theta->tau, theta->p, theta->tauv)!=0){
		fprintf(stderr,"Error in randomGen_init: the generated theta is\n");
		fprintf_THETA(theta,stderr,1);
		}

	return theta;
}


static void update_bic_solution(THETA *ptheta_bic, THETA theta_em)
{	THETA theta_bic = NULL;
	int flag1=0, flag2=0;
	theta_bic = *ptheta_bic;

	if(theta_em==NULL) return;
	
	if(theta_bic==NULL){
		theta_bic = duplicate_THETA(theta_em);
		}else{
		if(theta_em->pi<=1.0&&theta_em->pi>=0.0) flag2 = 1;
		if(theta_bic->pi<=1.0&&theta_bic->pi>=0.0) flag1 = 1;

		if((flag1==flag2&&theta_em->bic < theta_bic->bic) || (flag1==0&&flag2==1)){
			destroy_THETA(theta_bic); theta_bic = NULL;
			theta_bic = duplicate_THETA(theta_em);
			}
		}
	*ptheta_bic = theta_bic;

	return;
}

static void update_bic_solution_vec(THETA_VEC vec, THETA theta_em)
{	int flag1 = 0, flag2 = 0;


	if(theta_em==NULL) return;
	if(theta_em->ploidy > vec->lst || vec->theta[theta_em->ploidy]==NULL){
		set_value_THETA_VEC(vec, theta_em, theta_em->ploidy);
		} else{
		if(theta_em->pi<=1.0&&theta_em->pi>=0.0){ /*purity of this solution is within [0,1]*/
			flag2 = 1;
			}
		if(vec->theta[theta_em->ploidy]->pi<=1.0&&vec->theta[theta_em->ploidy]->pi>=0.0){ /*purity of this solution is within [0,1]*/
			flag1 = 1;
			}
		if((flag1==flag2&&theta_em->bic < vec->theta[theta_em->ploidy]->bic) || (flag1==0&&flag2==1)){
			set_value_THETA_VEC(vec, theta_em, theta_em->ploidy);
			}
		}
}


static THETA EM_auto(double *x, int B, int max_it, double epsilon, double max_ploidy, double an, double sx, THETA_VEC *vct_in)
{	THETA theta_init=NULL, theta_em=NULL, theta_bic=NULL, theta_tmp=NULL, theta_tmp1=NULL;
	double omega=0.1, tau0=0.1, sd0 = 0.1, D = 2.0;
	int alpha_fix = 0, flag=0, i;
	const double EPS = 1e-10;
	THETA_VEC vec = NULL;
	
	vec = *vct_in;
	if(vec==NULL) vec = new_THETA_VEC(5);


	flag = 0;
	for(omega=0.0;omega <= 1.0 + EPS; omega += 0.1){
		for(sd0= 0.1;sd0 <= 0.4 + EPS; sd0+=0.1){
			tau0 = sd0*sd0;
			for(D=2.0; D <= max_ploidy + EPS; D += 1.0){
				theta_init = create_init(x, B, omega, D, tau0,1, an, sx);
				theta_em =  new_THETA(theta_init->n, theta_init->EV, an, sx);
				if(theta_init->alpha==1.0) alpha_fix = 1; else alpha_fix = 0;
				EM_init(x,B,theta_init,theta_em, max_it, epsilon, alpha_fix);
				update_bic_solution_vec(vec, theta_em);
				destroy_THETA(theta_init);theta_init = duplicate_THETA(theta_em); theta_init->EV=0; /*prepare initial value for unequal variance case*/
				reduceModel(theta_em,x,B,max_it,epsilon); /*try to reduce the model*/
				update_bic_solution(&theta_bic,theta_em);
				update_bic_solution_vec(vec, theta_em);
				destroy_THETA(theta_em); theta_em = NULL;

				/*now unequl variance case*/
				theta_em = new_THETA(theta_init->n, theta_init->EV, an, sx);
				EM_init(x,B,theta_init,theta_em, max_it, epsilon, alpha_fix);
				update_bic_solution_vec(vec, theta_em);
				reduceModel(theta_em,x,B,max_it,epsilon); /*try to reduce the model*/
				update_bic_solution(&theta_bic,theta_em);
				update_bic_solution_vec(vec, theta_em);				

				destroy_THETA(theta_em); theta_em = NULL;
				destroy_THETA(theta_init); theta_init = NULL;
				}
			}
		}

	/*now try to reduce the models*/
	/*for(i=0;i<vec->size;i++){
		theta_tmp = duplicate_THETA(vec->theta[i]);
		if(theta_tmp!=NULL){
			reduceModel(theta_tmp,x,B,max_it,epsilon);
			update_bic_solution_vec(vec, theta_tmp);
			update_bic_solution(&theta_bic,theta_tmp);
			destroy_THETA(theta_tmp); theta_tmp = NULL;
			}
		}*/
	for(i=0;i<vec->size;i++){
		theta_tmp = duplicate_THETA(vec->theta[i]);
		if(theta_tmp!=NULL){
			theta_tmp1 = forceReduce(theta_tmp,x,B,max_it,epsilon);
			update_bic_solution_vec(vec, theta_tmp1);
			update_bic_solution(&theta_bic,theta_tmp1);
			destroy_THETA(theta_tmp); theta_tmp = NULL;
			destroy_THETA(theta_tmp1); theta_tmp1 = NULL;
			}
		}


	*vct_in = vec;
	return theta_bic;
}


static THETA EM_auto_new(double *x, int B, int max_it, double epsilon, int maxComp_in, double an, double sx, THETA_VEC *vct_in)
{	THETA theta_init=NULL, theta_em=NULL, theta_bic=NULL, theta_tmp=NULL, theta_tmp1 = NULL;
	int maxComp,  nc, alpha_fix, i;
	THETA_VEC vec=NULL;

	maxComp = maxComp_in;
	if(maxComp<4) maxComp=4;
	vec = *vct_in;
	if(vec==NULL) vec = new_THETA_VEC(maxComp+1);

	for(nc = 1; nc <= maxComp; nc++){
		theta_init = create_init_new(x, B, nc, 1, an, sx);
		theta_em = new_THETA(theta_init->n, 1, an, sx);
		if(theta_init->alpha==1.0) alpha_fix = 1; else alpha_fix = 0;
		EM_init(x,B,theta_init,theta_em, max_it, epsilon, alpha_fix);
		update_bic_solution_vec(vec, theta_em);

		/*try to simplify the model*/
		reduceModel(theta_em,x,B,max_it,epsilon);
		update_bic_solution(&theta_bic,theta_em);
		update_bic_solution_vec(vec, theta_em);

		
		/*unequal vairance case*/
		destroy_THETA(theta_init); theta_init = NULL;
		theta_init = new_THETA(theta_em->n,0, theta_em->an, theta_em->sx);
		theta_em->EV=0;
		copy_THETA(theta_em, theta_init);
		if(theta_init->alpha==1.0) alpha_fix = 1; else alpha_fix = 0;
		EM_init(x,B,theta_init,theta_em, max_it, epsilon, alpha_fix);
		update_bic_solution_vec(vec, theta_em);

		/*simplify the model*/
		reduceModel(theta_em,x,B,max_it,epsilon);
		update_bic_solution(&theta_bic,theta_em);
		update_bic_solution_vec(vec, theta_em);
		
		destroy_THETA(theta_em); theta_em = NULL;
		destroy_THETA(theta_init); theta_init = NULL;		
		}

	/*now try to reduce the models*/
	/*for(i=0;i<vec->size;i++){
		theta_tmp = duplicate_THETA(vec->theta[i]);
		if(theta_tmp!=NULL){
			reduceModel(theta_tmp,x,B,max_it,epsilon);
			update_bic_solution_vec(vec, theta_tmp);
			update_bic_solution(&theta_bic,theta_tmp);
			destroy_THETA(theta_tmp); theta_tmp = NULL;
			}
		}*/
	for(i=0;i<vec->size;i++){
		theta_tmp = duplicate_THETA(vec->theta[i]);
		if(theta_tmp!=NULL){
			theta_tmp1 = forceReduce(theta_tmp,x,B,max_it,epsilon);
			update_bic_solution_vec(vec, theta_tmp1);
			update_bic_solution(&theta_bic,theta_tmp1);
			destroy_THETA(theta_tmp); theta_tmp = NULL;
			destroy_THETA(theta_tmp1); theta_tmp1 = NULL;
			}
		}


	*vct_in = vec;
	return theta_bic;
}


static THETA EM_auto_randomstart(double *x, int B, int max_it, double epsilon, int maxComp_in, double an, double sx, int nRS, THETA_VEC *vct_in)
{	THETA theta_init=NULL, theta_em=NULL, theta_bic=NULL, theta_tmp=NULL, theta_tmp1=NULL;
	int maxComp,  nc, alpha_fix, EV,Nini,i;
	THETA_VEC vec=NULL;

	maxComp = maxComp_in;
	if(maxComp<4) maxComp=4;
	Nini = nRS;
	//if(Nini<1) Nini = 5;
	if(Nini<1) return NULL;
	vec = *vct_in;
	if(vec==NULL) vec = new_THETA_VEC(maxComp+1);

	for(nc = 1; nc <= maxComp; nc++){
		for(EV=0;EV<=1;EV++){
			for(i=0;i<Nini;i++){
				theta_init = randomGen_init(nc, EV, an, sx);
				theta_em = new_THETA(theta_init->n, 1, an, sx);
				if(theta_init->alpha==1.0) alpha_fix = 1; else alpha_fix = 0;
				EM_init(x,B,theta_init,theta_em, max_it, epsilon, alpha_fix);
				update_bic_solution_vec(vec, theta_em);
				//reduceModel(theta_em,x,B,max_it,epsilon);
				update_bic_solution(&theta_bic,theta_em);
				//update_bic_solution_vec(vec, theta_em);

				destroy_THETA(theta_em); theta_em = NULL;
				destroy_THETA(theta_init); theta_init = NULL;
				}
			}
		}


	/*now try to reduce the models*/
	/*for(i=0;i<vec->size;i++){
		theta_tmp = duplicate_THETA(vec->theta[i]);
		if(theta_tmp!=NULL){
			reduceModel(theta_tmp,x,B,max_it,epsilon);
			update_bic_solution_vec(vec, theta_tmp);
			update_bic_solution(&theta_bic,theta_tmp);
			destroy_THETA(theta_tmp); theta_tmp = NULL;
			}
		}*/

	for(i=0;i<vec->size;i++){
		theta_tmp = duplicate_THETA(vec->theta[i]);
		if(theta_tmp!=NULL){
			theta_tmp1 = forceReduce(theta_tmp,x,B,max_it,epsilon);
			update_bic_solution_vec(vec, theta_tmp1);
			update_bic_solution(&theta_bic,theta_tmp1);
			destroy_THETA(theta_tmp); theta_tmp = NULL;
			destroy_THETA(theta_tmp1); theta_tmp1 = NULL;
			}
		}


	*vct_in = vec;

	return theta_bic;
}


static double *filterOutlier(double *copyratio, int n, int *n_remain, THETA theta, double min_density, double min_comp_size)
{	int i, len, k;
	double *x=NULL;

	x = (double *) malloc(sizeof(double)*(n+10));
	if(x==NULL){
		fprintf(stderr,"Error in filterOutlier: memory allocation failed\n");
		}
	len = 0;
	for(i=0;i<n;i++){
		k = predict_copynumber(theta,copyratio[i],min_density);
		if(k!=-1 && theta->p[k] > min_comp_size){
			x[len] = copyratio[i];
			len++;
			}
		}
	*n_remain = len;

	return x;
}
