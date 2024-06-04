#ifndef purityEM_HX_XI
#define purityEM_HX_XI


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include "gamma.h"

typedef struct THETA_t{
	int n; /*maximum copy number*/
	double *p;
	int ploidy; /*the ploidy number; corresponding to the largest p[i]*/
	double alpha, xi, tau;
	double omega, D, pi;
	double *mu; /*auxilary variable: take values of 0, 1/2, 1, 3/2, 2, n/2*/
	double *mu_new; /*mu_new[i] = mu[i] + alpha*(mu[2]-mu[i]) + xi*/
	double *tauv; /*when tau is not a constant*/

	int *index_nonzero; /*index of the nonzero p_i*/
	int len_index,size_index; /*length of the index*/
	double num_param;

	int EV; /*whether equal variance (1) or not (0) */

	double bic, lk; /*update_THETA will not update these two values*/

	double an, sx; /*penalty coefficients for tiny values of tauv;*/

	double *psd, *pvalue; /*estimated standard deviation and pvalue of p*/
	int WF; /*the pvalue and psd were assigned (1) or not (0)*/
	} *THETA;
THETA new_THETA(int n, int EV, double an, double sx);
void copy_THETA(THETA theta1, THETA theta2); /*copy theta1 to theta2*/
THETA duplicate_THETA(THETA theta);
void destroy_THETA(THETA theta1);
int update_THETA(THETA theta, double alpha, double xi, double tau, double *p, double *tauv); 
/*update the values of theta as the given values, the length of p must be theta->n+1, though this function will NOT check; return nonzero if error*/
void fprintf_THETA(THETA theta, FILE *out, int header);
THETA del_component(THETA theta, int j); /*delete the jth component of theta and return a new one with the componet deleted*/
THETA setzero_component(THETA theta, int j); /*set the p[j] as 0.0 (note this is different from del_component)*/


double dnorm(double x, double mu, double sigma2); /*density function of normal distribution*/
double dnormMix(double x, double *mu, double *sigma2, double *p, int n); /*density of normal mixture model; this function will not check if sum(p)==1 and p>=0*/
double loglk(double *x, int B, THETA theta); /*calculate the log likelihood of the normal miture model*/
double regloglk(double *x, int B, THETA theta); /*calculate the regularized log likelihood of the normal miture model (with penalty of tiny tauvs)*/
double fisher_pi(double *x, int B, THETA theta, int i); /*calcuate partial^2 l/ partial p_i^2*/
void compPloidy(THETA theta); /*compute ploidy number, D, omega, pi, and num_param in theta*/

void Erho(double *rho, double *x, int B, int n, THETA thetat);
/*The E step of the EM algoirhtm, rho should be a B by (n+1) maxtrix
  thetat is the parameter estimate from the last M-step
  x is the copy ratio vector (length B)
  n should equal to thetat->n (though this function will not check for this)
*/
int Mstep(double *rho, double *x, int B, int n, THETA thetat, THETA thetat1, int alpha_fix); 
/*the M-step, the updated parameter will be stored at theta1
  return non-zero if error
  alpha_fix: alpha should be fixed in the iterations
 */

double diff_theta(THETA theta1, THETA theta2);


int EM_init(double *x, int B, THETA theta_ini, THETA theta_em, int max_it, double epsilon, int alpha_fix);
/* return nonzero if error
   theta_ini: the initial value of the EM algorithm
   thta_em: the resulting estimate from EM
*/
double peakHeight(THETA theta, int j); /*calculate the height of the ith compoent of theta*/
int which_min_Peak(THETA theta,int left); /*calculate the minimum peak*/
int assignTHETA_Pvalue(double *x, int B, THETA theta); /*return nonzero if error*/

THETA EM_purity(double *x, int B, int max_it, double epsilon); 
/* Assuming the ploidy number is 2,  
   this function automatically runs EM and estimate the purity starting from a series of initial values and return the estimate with minimum BIC
  here I assume x is sorted nondecreasingly and all of x are nonnegative
*/

int predict_copynumber(THETA theta, double x, double min_density);
/*predict the copy number of a given copy ratio x and a model theta
  return -1 if x does not belong to any component of the model (i.e. for each component, the density at x is smaller than min_density)
*/


void reduceModel(THETA theta, double *x, int B, int max_it, double epsilon); /*given a model, try to reduce the model complexity (if possible)*/
THETA forceReduce(THETA theta, double *x, int B, int max_it, double epsilon);


typedef struct THETA_VEC_t {
	THETA *theta;
	int lst,size; /*lst is the last non-NULL element in theta*/
	} *THETA_VEC;

THETA_VEC new_THETA_VEC(int n);
void destroy_THETA_VEC(THETA_VEC vc);
void extend_THETA_VEC(THETA_VEC theta_v, int n); /*make the vector larger*/
void set_value_THETA_VEC(THETA_VEC theta_v, THETA theta, int i); /*set the ith element of theta_v as theta*/

#endif
