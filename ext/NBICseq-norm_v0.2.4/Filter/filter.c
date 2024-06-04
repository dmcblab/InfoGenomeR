/*
Written by Ruibin Xi
*/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <getopt.h>
#include <math.h>
#include <unistd.h>
#include "../lib/read.h"

typedef struct Interval_t{
        int start, end;
        } Interval;

typedef struct IntervalSet_t{
	Interval *I;
	int len, size;
	} *IntervalSet;

IntervalSet newIntervalSet();
void append_Interal(IntervalSet I, int start, int end);
void fprintf_InvervalSet(FILE *out,IntervalSet I);
int cmpInterval(const void *a , const void *b); /*if a and b have any overlap, they will be set as equal*/

static int is_number(char c);
/*test whether a character c is '0','1',..,'9' or '-' (minus sign)
 *  *if c is '0','1',..,'9' return 1
 *   *if c is '-' return -1
 *    *otherwise return 0
 *     * */

static int is_numeric(char *ll);
/*test if a character array starts with numeric values
 *return 1 if true, 0 otherwise
 * */


int readOneCNV(FILE *in, char *chrom, IntervalSet Iset); 
/*chrom is the given chromosome name
  return nonzero if error
*/ 

int readCNVonChrom(FILE *in, char *chrom, IntervalSet Iset);



int main(int argc, char ** argv)
{	FILE *cnvRegion=NULL, *readCntFile=NULL;
	char *chrom=NULL, *ltxt;
	IntervalSet Iset=NULL;
	Interval Itmp;
	int pos;
	//int nrow, ncol;
	//double *dta=NULL;

	Iset = newIntervalSet();
	if(argc!=4){
		fprintf(stderr,"Usage: %s <cnvRegionFile> <readCountFile> <ChromName>\n", argv[0]);
		exit(1);
		}

	cnvRegion = fopen(argv[1],"r");
	if(cnvRegion==NULL){
		fprintf(stderr,"No such file or directory: %s\n", argv[1]);
		exit(1);
		}
	readCntFile = fopen(argv[2],"r");
	if(readCntFile==NULL){
		fprintf(stderr,"No such file or directory: %s\n", argv[2]);
		exit(1);
		}
	chrom = strdup(argv[3]);

	readCNVonChrom(cnvRegion,chrom, Iset);
	qsort(Iset->I,Iset->len, sizeof(Interval), cmpInterval);
	//fprintf_InvervalSet(stdout,Iset);

	while(!feof(readCntFile)){
		//dta = read_table(readCntFile, &nrow, &ncol, 1, 0);
		ltxt = lr_read(readCntFile);
		pos = atoi(ltxt); /*this should be the first item in the line, i.e. the position on the chromosomes*/
		//flag = 0; /*not a CNV region*/
		//if(nrow==1&&ncol==2&&dta!=NULL){
		if(ltxt!=NULL){
			Itmp.start = Itmp.end = pos; //(int) dta[0];
			if(bsearch(&Itmp, Iset->I, Iset->len, sizeof(Interval),cmpInterval)==NULL){ /*this position is not in the CNV region*/
				fprintf(stdout,"%s",ltxt);
				}
			//if(flag==1){fprintf(stdout,"%g\t%g\n",dta[0],dta[1]); getchar();}
			//if(flag==1){fprintf(stdout,"%s",ltxt);}
			}

		/*if(dta!=NULL){
			free(dta);
			dta=NULL;
			}*/
		}


	return 0;
}

IntervalSet newIntervalSet()
{	IntervalSet Iset = malloc(sizeof(struct IntervalSet_t));
	Iset->size = 10;
	Iset->len = 0;
	Iset->I = (Interval *) malloc(sizeof(Interval)*Iset->size);
	if(Iset==NULL){
		fprintf(stderr,"Error in newIntervalSet: memory allocation failed\n");
		exit(1);
		}	

	return Iset;
}


static int is_number(char c)
        { if(c<='9'&&c>='0') return 1;
          else if(c=='-') return -1;
          else return 0;
        }

static int is_numeric(char *ll)
        { int len,flag;
          len = strlen(ll);

          if(len==0) return 0;

          flag = is_number(ll[0]);
          if(flag==0) return 0;
          else if(flag==-1)
                { if(len==1||is_number(ll[1])!=1) return 0;
                  else return 1;
                }
          else return 1;
        }


void append_Interal(IntervalSet Iset, int start, int end)
{
	if(Iset->len +2 > Iset->size){
		Iset->size += 100;
		Iset->I = (Interval *) realloc(Iset->I, sizeof(Interval)*Iset->size);
		if(Iset->I==NULL){
			fprintf(stderr,"Error in append_Interal: memory allocation failed\n");
			exit(1);
			}
	}

	Iset->I[Iset->len].start = start;
	Iset->I[Iset->len].end = end;
	Iset->len ++;
	
	return;
}

void fprintf_InvervalSet(FILE *out,IntervalSet Iset)
{	int i;
	for(i=0;i<Iset->len;i++){
		fprintf(out,"%d\t%d\n",Iset->I[i].start,Iset->I[i].end);
		}
	return;
}


int readOneCNV(FILE *in, char *chrom, IntervalSet Iset)
{	char *ll=NULL, *sub_str=NULL;
	int start, end;
	
	ll = lr_read(in);
	
	if(strlen(ll)==0) return 0;
	sub_str = strtok(ll," \t"); /*The first column should be the chromosome name*/
	//fprintf(stderr,"substr=%s\n",sub_str);
	//fprintf(stderr,"chrom=%s\n",chrom);
	if(strcmp(chrom,sub_str)!=0) return 0; /*not on the required chromosome*/
	//else{fprintf(stderr,"Yes\n");}
	/*2nd column*/
	sub_str=strtok(NULL," \t");
	if(sub_str!=NULL && is_numeric(sub_str)){
		start = atoi(sub_str);
		}else{
		return 0;
		}
	/*3rd column*/
	sub_str = strtok(NULL," \t");
        if(sub_str!=NULL && is_numeric(sub_str)){
                end = atoi(sub_str);
                }else{
                return 0;
                }

	append_Interal(Iset, start, end);

	return 0;
}


int readCNVonChrom(FILE *in, char *chrom, IntervalSet Iset)
{	int i=0;
	while(!feof(in)){
		readOneCNV(in, chrom,Iset);
		i++;
		}

	//fprintf(stderr,"Scanned %d lines; %d regions on chromosome %s\n", i, Iset->len, chrom);
	return 0;

}



int cmpInterval(const void *a , const void *b)
        { Interval *I1, *I2;
          I1 = (Interval *) a;
          I2 = (Interval *) b;

          if(I1->end < I2->start){return -1;}
          else if(I1->start > I2->end){return 1;}
          else return 0;
        }

