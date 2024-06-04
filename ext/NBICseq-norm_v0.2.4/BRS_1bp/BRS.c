#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <time.h>
#include <getopt.h>
#include "../lib//read.h"
#include "pos_cnt_lst.h"

#define GAP_ERROR 2

static int min_dist_to_gap = 1000000;
static int tt_read_mp = 0; /*total reads processed on the mappable positions*/
static int tt_mappable = 0;

typedef struct SRM_binning_t{
	char *tumor_file,*mapfile,*output, *gap_file;
	int win_size; /*the window size for removing the singular reads*/
	double quantile, multple;
	FILE *file_sum; 
	} SRM_binning; /*the argument passed to the program*/

typedef struct RDcountBase_t{
	int pos, count;
	} RDcountBase; /*a data structure that used to store the read count at certain genomic position*/

typedef struct Interval_t{
        int start, end;
        } Interval;


static int cmpdouble(const void *a, const void *b);

static int cmpRDcountBase(const void *a, const void *b);

static int *aggregate(double *reads, int n_reads, int *num_pos);
/* aggregate the reads located at the same position together.
 * Returned value is a two column matrix looking like
 * 
 * <position>	<number of read at this position>
 *
 * The argument num_pos will record the size of this matrix (number of rows)
 *
 * reads: the read positions which should be ordered nondecreasingly
 * n_reads: total number of reads (length of the vector reads)
 * */


static int *sort_rms_binning(double *tumor, int n_tmor, int *num_tum_1bp,int w, double quantile, double multple);
/* sort the reads, remove the singular reads and bin the data to 1bp bin.
   arguments:
	num_tum_1bp: the total number of bins obtained.
	w: the window size used to determine the singular genomic positions.
	quantile: quantile (e.g. 0.95) used to determine the singular genomic positions.
	multple: if a genomic position has more than multple * quantile, it will be identified as an outlier.
 */
RDcountBase *create_RDconntBase_From1bpbin(int *bin, int nbins); /*bin: 1 bp bin, nbins: number of 1 bp bins*/


Interval * read_gapfile(char *gapfile, int *num_interval);
int cmpInterval(const void *a , const void *b);

int distanc_gap(int pos, Interval *gaps, int ngaps); /*distance in kb*/

static void bin_fprint(int *bin,int nbins,FILE *output);

static void RDcountBase_fprint(RDcountBase *RDcount, int nbins,double *map, int nrow_map,FILE *output, Interval *gaps, int ngaps);
/* RDcount should have been sorted according to pos
*/

static void explain_command(char **argv);
static SRM_binning option_assign(int argc, char **argv);



int main(int argc, char **argv)
{
	int n_tmor,ncol,nbins,*bin, nrow_map,ncol_map;
	FILE *in_tmor,*output, *inmap;
	double  *tumor=NULL, *map=NULL;
	SRM_binning args;
	RDcountBase *bin_RD=NULL; 
	Interval *gaps=NULL;
	int ngap=0;
	
	args = option_assign(argc,argv);

	if(args.gap_file!=NULL){
		gaps = read_gapfile(args.gap_file,&ngap);
		}
	
	in_tmor = fopen(args.tumor_file,"r");
	if(in_tmor == NULL)
		{ fprintf(stderr,"fopen %s: %s\n",args.tumor_file ,strerror(errno));
		  exit(1);
		}

	if(args.mapfile!=NULL){
		inmap = fopen(args.mapfile,"r");
		if(inmap==NULL){
			fprintf(stderr,"fopen %s: %s\n",args.mapfile,strerror(errno));
			exit(1);
			}
		map = read_table(inmap,&nrow_map,&ncol_map,-1,0);
		if(ncol_map!=2){
			fprintf(stderr,"Error, mappability file should be a 2 column data file\n");
			exit(1);
			}
		}	


	tumor = read_table(in_tmor,&n_tmor,&ncol,-1,0);
	if(ncol!=1) {fprintf(stderr,"Error: seq file has multiple columns.\n"); exit(1);}
	fprintf(stderr,"%d reads loaded\n",n_tmor);

	if(tumor==NULL) {fprintf(stderr,"No short read positions loaded.\n"); exit(0);} 

	bin = sort_rms_binning(tumor,n_tmor,&nbins,args.win_size,args.quantile,args.multple);
	/*sort, remove singular positions and bin*/


	if(args.mapfile!=NULL){
		bin_RD = create_RDconntBase_From1bpbin(bin,nbins); /*since bin has been sorted, no need to sort bin_RD again*/
		free(bin);
		bin = NULL;
		}

	if(args.output!=NULL){
        	output = fopen(args.output,"w");
        	if(output==NULL) {fprintf(stderr,"fopen %s: cannot create the file.\n",args.output); exit(1);}

		if(args.mapfile==NULL) {bin_fprint(bin,nbins,output); free(bin);bin=NULL;}
		else { 
			RDcountBase_fprint(bin_RD, nbins, map, nrow_map, output, gaps, ngap); free(bin_RD);
			if(args.file_sum!=NULL){
				fprintf(args.file_sum,"Total_Reads_On_Mappable_Positions\tNum_Mappable_Positions\n");
				fprintf(args.file_sum,"%d\t%d\n",tt_read_mp,tt_mappable);
				}
			}
		}
	else {
		if(args.mapfile==NULL){
			bin_fprint(bin,nbins,stdout);
			free(bin); bin = NULL;
			}else{
			RDcountBase_fprint(bin_RD, nbins, map, nrow_map, stdout, gaps, ngap); 
			free(bin_RD); bin_RD = NULL;
			if(args.file_sum!=NULL){
				fprintf(args.file_sum,"Total_Reads_On_Mappable_Positions\tNum_Mappable_Positions\n");
				fprintf(args.file_sum,"%d\t%d\n",tt_read_mp,tt_mappable);
				}
			}
		}
	return 0;	

}



static void explain_command(char **argv)
	{ fprintf(stderr,"Usage: %s [options] <seq file>\n",argv[0]);
	  fprintf(stderr,"<seq file>: a one column file only containing the positions of short reads.\n");
	  fprintf(stderr,"Options:\n");
	  fprintf(stderr,"-h: print this message.\n");
	  fprintf(stderr,"-q <float>: the quantile used for identification of the singular genomic positions; default is 0.95\n");
	  fprintf(stderr,"-w <int>: the window size for calculating the quantiles; default is 200\n");
	  fprintf(stderr,"-o <string>: the output bin file; if unspecified, print to the stdout.\n");
	  fprintf(stderr,"-m <string>: the mappability file; If provided, filter by mappability and the output will contain uniquely mappable positions with no reads mapped\n");
	  fprintf(stderr,"-g <string>: the gap file (two column data file with 1st column the start positions of the gaps and the 2nd column the end positions of the gaps)\n");
	  fprintf(stderr,"-s <string>: report the summary statistics to the file <string>\n");
	  fprintf(stderr,"--multiplicity <float>: If a genomic position has more than multiplicity*quantile number of reads,\n");
	  fprintf(stderr,"                        it will be viewed as an outlier\n");
	  fprintf(stderr,"                        and the number of reads at this position will be set as multiplicity*quantile;\n");
	  fprintf(stderr,"                        default is 5.0\n");
	  return;
	}


static SRM_binning option_assign(int argc, char **argv)
	{ SRM_binning args;
	  char *file_sum = NULL;
	  int c,option_index;
	  static struct option long_options[] =
			{ {"multiplicity", required_argument ,0},
			  {0,0,0,0} 
			};

	  args.multple = 5.0;
	  args.quantile = 0.95;
	  args.win_size = 200;
	  args.output=NULL;
	  args.gap_file = NULL;
	  args.mapfile = NULL;
	  args.file_sum = NULL;

	  while((c=getopt_long(argc,argv,"g:m:q:b:w:o:s:h",long_options,&option_index))!=-1){
		switch(c){
			case 0:
				args.multple = atoi(optarg);
				if(args.multple<1.0) { fprintf(stderr,"Error,multiplicity must be larger than or equal to 1.0\n"); exit(1);}
				break;
			case 'q':
				args.quantile = atof(optarg);
				if(args.quantile<=0.0 || args.quantile>1.0)  { fprintf(stderr,"Quantile should be between 0 and 1.\n"); exit(1);}
				break;
			case 'w':
				args.win_size = atoi(optarg);
				if(args.win_size<1) {fprintf(stderr,"The window size must be positive.\n");exit(1);}
				break;
			case 'g':
				args.gap_file = strdup(optarg);
				break;
			case 's':
				file_sum = strdup(optarg);
				break;
			case 'o':
				args.output = strdup(optarg);
				break;
			case 'm':
				args.mapfile = strdup(optarg);
				break;
			case 'h': 
				explain_command(argv);
				exit(0);
			case '?': /* getopt_long already printed an error message. */
				exit(1);
				break;
			default:
				abort (); 
			}
		}

        if (argc - optind!=1){
           	explain_command(argv); exit(1);
         	}

	args.tumor_file = strdup(argv[optind]);

	if(file_sum!=NULL){
		args.file_sum = fopen(file_sum,"w");
		if(args.file_sum==NULL){
			fprintf(stderr,"fopen %s: %s\n",file_sum,strerror(errno));
			exit(1);
			}
		free(file_sum);
		file_sum=NULL;
		}

	return args;

	}


static int cmpdouble(const void *a, const void *b)
	{ double tmp1,tmp2;
	  tmp1 = *((const double *) a);
	  tmp2 = *((const double *) b);
	  if(tmp1<tmp2) return -1;
	  else if(tmp1>tmp2) return 1;
	  else return 0;
	}

static int cmpRDcountBase(const void *a, const void *b)
	{RDcountBase tmp1, tmp2;
	tmp1 = *((RDcountBase *)a);
	tmp2 = *((RDcountBase *)b);
	if(tmp1.pos<tmp2.pos) return -1;
	else if(tmp1.pos>tmp2.pos) return 1;
	else return 0;	
	}


static int *aggregate(double *reads, int n_reads, int *num_pos)
	{ /*The first column of read_dist is the position of the read*/
	  /*The second column of read_list is the number of reads at this position.*/
	  /*If certain position does not have any reads, this position will be ignored*/
	  int *reads_dist, pos_min, pos_max,size;
	  int i,j;
	
	  pos_min = (int) reads[0];
	  pos_max = (int) reads[n_reads-1];

	  size = (n_reads<pos_max-pos_min+1)? n_reads : pos_max-pos_min+1;
	  reads_dist = (int *)  malloc(sizeof(int)*(size*2+10));

	  j = 0;
	  i =0;
	  reads_dist[2*j] = (int) reads[i]; /*The first column of read_dist is the position of the read*/
	  reads_dist[2*j+1] = 1; /*The second column of read_list is the number of reads at this position.*/

	  for(i=1;i<n_reads;i++)
		{ if((int) reads[i]==reads_dist[2*j]) 
			{ reads_dist[2*j+1]++;} /*one more read at position reads_list[2*j]*/
		  else /*get to a new position*/
			{ j++;
			  reads_dist[2*j] = (int) reads[i];
			  reads_dist[2*j+1] = 1;
			}
		}

	 *num_pos = j+1;
	
	 return reads_dist;

	}


RDcountBase *create_RDconntBase_From1bpbin(int *bin, int nbins)
	{int i;
	RDcountBase *rdcnt=NULL;
	if(nbins<=0){
		return NULL;
		}
	rdcnt = (RDcountBase *) malloc(sizeof(RDcountBase)*(nbins+1));
	if(rdcnt==NULL) {fprintf(stderr,"Error in create_RDconntBase_From1bpbin: memory allocation failed\n");exit(1);}

	for(i=0;i<nbins;i++){
		rdcnt[i].pos = bin[2*i];
		rdcnt[i].count = bin[2*i+1];
		}
	return rdcnt;
	}


static int *sort_rms_binning(double *tumor, int n_tmor, int *num_tum_1bp,int w, double quantile, double multple)
	{ int *tumor_1bpbin;

	  qsort(tumor,n_tmor,sizeof(double),cmpdouble);
          fprintf(stderr,"sorted %d reads\n",n_tmor);
 
	  /*aggregate the reads; */
	  tumor_1bpbin = aggregate(tumor,n_tmor,num_tum_1bp);
	  free(tumor); /*If invoke in R, it's maybe better not free tumor here*/

	  /*remove the singular points*/
	  singularity_rm(tumor_1bpbin, *num_tum_1bp, w, quantile, multple);


	 return tumor_1bpbin;
	}


static void bin_fprint(int *bin,int nbins,FILE *output)
{	int i;

        if(bin==NULL||nbins<=0) {fprintf(stderr,"bin_print: empty array.\n"); return;}
        for(i=0;i<nbins;i++)
		{ 
	 	fprintf(output,"%d\t",bin[2*i]);
		fprintf(output,"%d\n",bin[2*i+1]);
                }
        return;

}




static void RDcountBase_fprint(RDcountBase *RDcount, int nbins,double *map, int nrow_map,FILE *output, Interval *gaps, int ngaps)
{	RDcountBase tmp1, *search_rslt=NULL;
	int i,pos;
	int start, end;
	int dist;
	
	tmp1.count = 0;
	for(i=0;i<nrow_map;i++){
		start = (int) map[2*i];
		end = (int) map[2*i+1];
		for(pos=start;pos<=end;pos++){
			tmp1.pos = pos;
			dist = distanc_gap(pos, gaps, ngaps);
			if(dist > min_dist_to_gap){
				tt_mappable++;
				search_rslt = bsearch(&tmp1,RDcount,nbins,sizeof(RDcountBase),cmpRDcountBase);
				if(search_rslt==NULL){
					fprintf(output,"%d\t%d\n",pos,0);
					}else{
					fprintf(output,"%d\t%d\n",search_rslt->pos,search_rslt->count);
					tt_read_mp += search_rslt->count;
					}
				search_rslt = NULL;
				}
			}	
		}	
	return;
}


Interval *read_gapfile(char *gapfile, int *num_interval)
        { Interval *Is=NULL;
          double *dta=NULL;
          int nrow,ncol,i;
          FILE *fin=NULL;

          fin = fopen(gapfile,"r");
          if(fin==NULL){
                fprintf(stderr,"No such file or directory: %s\n",gapfile);
                exit(1);
                }

          dta = read_table(fin,&nrow,&ncol,-1,0);
          if(nrow>0&&ncol!=2){
                fprintf(stderr,"The gap file must of 2 column\n");
                exit(GAP_ERROR); /*2 for gap file errors*/
                }
         if(nrow>0){
                Is = (Interval *) malloc(sizeof(Interval)*(nrow+1));
                if(Is==NULL){
                        fprintf(stderr,"memory allocation faile\n");
                        exit(1);
                        }
                *num_interval = nrow;
                }else{
                *num_interval = 0;
                }

        for(i=0;i<nrow;i++){
                Is[i].start = (int) dta[2*i];
                Is[i].end = (int) dta[2*i+1];
                if(Is[i].start>Is[i].end||Is[i].start<0){
                        fprintf(stderr,"Unexpected gap information in %s\nstart\tend\n%d\t%d\n",gapfile,Is[i].start,Is[i].end);
                        exit(GAP_ERROR); /*2 for gap file errors*/
                        }
                }
        if(dta!=NULL) {free(dta);dta=NULL;}
	if(nrow>0) {qsort(Is,nrow,sizeof(Interval), cmpInterval);}

        return Is;
        }


int cmpInterval(const void *a , const void *b)
        { Interval *I1, *I2;
          I1 = (Interval *) a;
          I2 = (Interval *) b;

          if(I1->start<I2->start){return -1;}
          else if(I1->start>I2->start){return 1;}
          else{
                if(I1->end<I2->end) return -1;
                else if(I1->end>I2->end) return 1;
                else return 0;
                }
        }

int distanc_gap(int pos, Interval *gaps, int ngaps)
        { int i;
          int min_dist=0, dist;

	  min_dist = min_dist_to_gap+10;
          if(ngaps>0 && gaps!=NULL ){
                for(i=0;i<ngaps;i++){
                        if(pos<gaps[i].start) dist =  gaps[i].start - pos;
                        else if(pos>gaps[i].end) dist = pos - gaps[i].end;
                        else dist = 0;

                        if(i==0) min_dist = dist;
                        min_dist = (min_dist<dist ? min_dist : dist);
                        }
                }

        return min_dist;
        }



