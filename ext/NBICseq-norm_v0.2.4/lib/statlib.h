/*This is the statistics library head file writing by Ruibn Xi Oc. 15 2009
 * Modified on Dec 16 2009
 * add a fucntion void db_shuffle(double *array, int n);
 * */

#ifndef STATLIB_XI
#define STATLIB_XI

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void seed_set(int seed); /*reset the seed of the random number generator rand_lp as seed.*/ 
double rand_lp();
/*Long period (> 2 * 10^18) random number generator of L'Ecuyer with Bays-Durham shuffle and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the end point values). Call with idum a negative integer to initialize; thereafter, do not alter idum between successive deviates in a sequence. RNMX should approximate the largest floating value that is less than 1
See Numerical Recipes in C++, page 286.
*/

void db_shuffle(double *array, int n);
void db_shuffle_int(int *array, int n);

double rgamma1(double alpha, double beta);
/* Generate Gamma distribution using the algorithm Ahrens and Dieter (1974)
 * See Monte Carlo Statistical Methods (C. Robert and G. Casella 2004) Page 22
 * alpha: the shape parameter
 * beta: the scale parameter
 * the density of gamma distribution is f(x;alpha,beta) = x^{alpha-1}\frac{e^{-x/beta}}{beta^{alpha}Gamma(alpha)}
 * return -1, if error and print an error message
 * */

int rDirichlet(double *alpha, double *x,int n);
/* Generate random variates from  Dirichlet distribution with parameter alpha,
 * both alpha and x shoulb be of length n
 * if error, return -1;
 * The gererated value will be stored in x;
 * */

#endif
