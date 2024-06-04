/*This is the statistics library written by Ruibin Xi Oct 15 2009*/


#include "statlib.h"
static int idum=-849345; /*the seed of the random generator rand_lp*/

void seed_set(int seed)
	{ idum = seed;}

double rand_lp()
	{ const int IM1 = 2147483563, IM2 = 2147483399;
	  const int IA1 = 40014, IA2 = 40692, IQ1 = 53668, IQ2 = 52774;
	  const int IR1 = 12211, IR2 = 3791, NTAB = 32, IMM1 = IM1-1;
	  const int NDIV = 1+IMM1/NTAB;
	  const double EPS = 3.0e-16, RNMX = 1.0- EPS, AM = 1.0/((double) IM1);
	  static int idum2 = 123456789, iy=0;
	  static int iv[32];
	  int j,k;
	  double temp;

	  //iv = (int *) malloc(sizeof(int)*(NTAB)); // initialize iv as a vector

	  if(idum<=0){ //Initialize
		  idum = (idum==0?1:-idum);
		  idum2 = idum;
		  for(j=NTAB+7;j>=0;j--){ // Load the shuffle talbe (after 8 warm ups)
			  k = idum/IQ1;
			  idum = IA1*(idum-k*IQ1)-k*IR1;
			  if(idum<0) idum += IM1;
			  if(j< NTAB) iv[j] = idum;
			}
		  iy = iv[0];
		}
	  k = idum/IQ1; // Start here when not initializing.
	  idum = IA1*(idum-k*IQ1)-k*IR1;
	  if(idum<0) idum +=IM1;
	  k = idum2/IQ2;
	  idum2=IA2*(idum2-k*IQ2)-k*IR2;
	  if(idum2<0) idum2 +=IM2;
	  j = iy/NDIV;
	  iy = iv[j]-idum2;
	  iv[j] = idum;
	  if(iy<1) iy += IMM1;
	  if((temp=AM*iy)>RNMX) return RNMX;
	  else return temp;
	}



void db_shuffle(double *array, int n)
	{ int i,k;
	  double tmp;
	  k = n;
	  while(k>1)
		{ 
		  i = (int) floor(k*rand_lp());
		  k--;
		  tmp = array[k];
		  array[k] = array[i];
		  array[i] = tmp;
		}
	  return;
	}

void db_shuffle_int(int *array, int n)
        { int i,k;
          int tmp;
          k = n;
          while(k>1)
                {
                  i = (int) floor(k*rand_lp());
                  k--;
                  tmp = array[k];
                  array[k] = array[i];
                  array[i] = tmp;
                }
          return;
        }


double rgamma1(double alpha, double beta)
	{ double U0, U1;
	  double x, y;
	  double e = 2.718281828459045; /*Euler's number*/
	  
	  if(alpha<=0||beta<=0) {printf("Error in rgamm1: invalid parameters.\n"); return -1;}

	  if(alpha==1) /*when alpha=1, gamma distribution is exponential distribution; simple generating method can be used*/
		{ U0 = rand_lp();
		  x = -log(1-U0);
		  return beta*x;
		}


	  do{
		  U0 = rand_lp();
	          U1 = rand_lp();
	  	  if(U0>e/(e+alpha))
			{ x = -log((alpha+e)*(1-U0)/alpha) + 1;
		  	  y = pow(x,alpha-1);
			}
	  	  else
			{ x = pow((alpha+e)*U0/e,1/alpha);
		  	  y = exp(-x);
			}
		}while(U1>=y);

	  return beta*x;
	}



int rDirichlet(double *alpha, double *x,int n)
	{ double sum;
	  int i;
	  
	  if(alpha==NULL||x==NULL) {printf("Error in rDirichlet.\n"); return -1;}
	  sum=0.0;
	  for(i=0;i<n;i++)
		{ x[i] = rgamma1(alpha[i],1.0);
		  if(x[i]<0) return -1;
		  sum = sum+x[i];
		}
	  for(i=0;i<n;i++)
		{ x[i] = x[i]/sum;}
	  return 1;
	}
