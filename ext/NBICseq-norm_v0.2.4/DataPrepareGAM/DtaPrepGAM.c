/* Written by ruibin Xi
   May 19, 2011
   Modified on Feb 4 2011
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include "../lib/read.h"
#include "../lib/statlib.h"
#include <getopt.h>
#include <math.h>

#define GAP_ERROR 2

static int no_header=0; /*if 1, do NOT print header*/

typedef struct DPargs_t{
        char *fa_file_name, *readcnt_file_name,*output_name,*chrom_name;// *gap_file;
        int fragment_size, readLen; /*the window size for removing the singular reads*/
	int extd_bp; /*how many basepairs to extend to upstream and downstream*/
	int nogapInRead; /*when use every nucleotide in [-extd_bp, readLen+extd_tp] instead of just the nucleotides on two ends*/
	int NumNuclotide; /*How many nucloetide considered in the glm model*/

        int estimate; /*if 1, try to estimate the expected count*/
	int bin_size;
	//int GC_power;
	int gc_bin; /*if 1, report the read count and gc in bins;*/
	//int gc_log; /*if 1 use log GC as predictor*/
	//double gc_mp; /*use gc + gc*gc^gc_mp instead of gc+gc*gc (in case of GC_power = 1)*/ 
	int negBinom; /*if 1 use negBinom, otherwise no*/
	int map_bin;

	int refine; /*if 1, refine the expected value with the given coefficient of GC*/
	double *rf_gc_coef; /*gc coeffients in the refinement step*/
	int size_rf_gc_coef;

	FILE *input_pois; /*the file contains the estimates of the poisson estimates*/
	FILE *fa_file, *rdcnt_file, *output;
	
	double inverse_theta; /*the inverse of theta parameter in the negative binomail model obtained from input_pois*/ 
        } DPargs; /*the argument passed to the program*/

typedef struct Interval_t{
	int start, end;
	} Interval;



static DPargs option_assign(int argc, char **argv);


char *read_fa(FILE *input,int *len,char **chrom); /*len will be the length of the chomosome*/


void print_header(DPargs args);
void print_oneline(DPargs args, char *dna, int dna_len, int pos,int read_count );

int count_nuc(char *dna, int length, int start, int end, char c);

static int update_gc(char *dna, int dna_len, int pos, int pos_pre, int fragment_size, int gc_pre);




double calculate_bias(DPargs args, char *dna, int dna_len, int pos, double *param_est, int nrow_param, int ncol_param, double *GC_bias);
/*return -1.0 if not calculable for the position*/


void fprint_gc_bin(FILE *output, char*dna, int dna_len, double *rdcnt, int nrow_rdcnt,  int bin_size, int map_bin, char *chrom_name); /*map_bin=1 use mappability bin*/

void fprintf_predicted(DPargs args, double *rdcnt, int nrow_rdcnt,  char *dna, int dna_len, double *param_est, int nrow_param, int ncol_param);


int main(int argc, char **argv)
{	DPargs args;
	char *dna=NULL, *chrom=NULL, *chrom_name=NULL; /*chrom_name is the chromosome name provided by the user*/
	int i,dna_len, pos,read_count;
	double *rdcnt=NULL;
	int ncol_rdcnt=0, nrow_rdcnt=0,nrow_param, ncol_param;
	double *param_est=NULL;

	

	args = option_assign(argc,argv);




	rdcnt = read_table(args.rdcnt_file,&nrow_rdcnt,&ncol_rdcnt,-1,0);
	if(ncol_rdcnt!=2){
		fprintf(stderr,"The ReadCountFile must be a two column file\n");
		exit(1);
		}


	dna = read_fa(args.fa_file,&dna_len,&chrom);
	if(chrom[strlen(chrom)-1]=='\n') chrom[strlen(chrom)-1] = '\0';
	fprintf(stderr,"Finished reading the fa file %s\n", args.fa_file_name);


	if(chrom!=NULL&&strlen(chrom)!=0&&args.chrom_name!=NULL&&strcmp(chrom,args.chrom_name)!=0) 
		fprintf(stderr,"Warning: the chromosome name given by the user (%s) and the chromosome name (%s) in %s are different\n",chrom_name,chrom,args.fa_file_name);


	if(args.estimate==1){
		param_est = read_table(args.input_pois,&nrow_param,&ncol_param,-1,2);
		if(args.negBinom==1){
			if(nrow_param<=1) {fprintf(stderr,"Incorrect format of parameter estimate file\n");exit(1);}
			if(param_est[(nrow_param-1)*ncol_param] >0.0) {args.inverse_theta = 1.0/param_est[(nrow_param-1)*ncol_param];}
			else {fprintf(stderr,"Warning,invalide estimate of the theta parameter = %g <= 0; Ignored\n",param_est[nrow_param*ncol_param]);}
			fprintf(stderr,"theta = %g\t 1/theta = %g\n",param_est[(nrow_param-1)*ncol_param], args.inverse_theta);
			nrow_param --;
			}
		if(nrow_param<=0) {fprintf(stderr,"Incorrect format of parameter estimate file\n");exit(1);}
		if(nrow_param != 3*args.NumNuclotide +  2*args.fragment_size + 2) { /*pay attention to here*/
			fprintf(stderr,"Error, there are %d covariats but %d parameter estimates\n",3*args.NumNuclotide + 2*args.fragment_size + 2, nrow_param); exit(1);
			}
		}

	if(args.estimate!=1){
		if(args.gc_bin==0){
			if(no_header!=1) print_header(args);
			for(i=0;i<nrow_rdcnt;i++){
				pos = (int) rdcnt[2*i];
				read_count = (int) rdcnt[2*i+1];
				if(pos<0||pos>dna_len){
					fprintf(stderr,"No such position (%d) in chrom %s\n",pos,chrom);
					exit(1);
					}
				if(read_count<0){
					fprintf(stderr,"Error, negative read count (%d), line %d in %s\n",read_count,i,argv[1]);
					}
				print_oneline(args,dna,dna_len,pos,read_count);
				}
			}else{
			fprint_gc_bin(args.output, dna,dna_len,rdcnt,nrow_rdcnt,args.bin_size, args.map_bin,args.chrom_name);
			}
		}else{
		fprintf(stderr,"Bin size = %d\n",args.bin_size);
		fprintf_predicted(args,rdcnt,nrow_rdcnt,dna,dna_len,param_est,nrow_param,ncol_param);
		}

	fclose(args.fa_file);

	return 0;
}





static void explain_command(char **argv)
        { fprintf(stderr,"usage: %s <fa file>\n",argv[0]);
	  fprintf(stderr,"The ouput will be printed on stdout\n");
	  fprintf(stderr,"Options:\n");
	  fprintf(stderr,"       -i <string>: a two column read count file(position and count); If unspecified, use <stdin> as input;\n");
	  fprintf(stderr,"                    Data should be ordered by position\n");
          fprintf(stderr,"       -s <int>: fragment size as estimated from paired end data or by cross correlation of reads on postive and negative strands; Default 300\n");
	  fprintf(stderr,"       -l <int>: read length; Default 50\n");
	  fprintf(stderr,"       -o <string>: if specified, use this as output\n");
	  fprintf(stderr,"       -e <string>: the file that stores the estimates of the glm model\n");
	  fprintf(stderr,"       -b <int>: the bin size; default 100; only valid if -e is specified\n");
          fprintf(stderr,"       -h: print this message.\n");
	  fprintf(stderr,"       --chrom <string>: the chromosome under consideration\n");
	  fprintf(stderr,"       --uds <int>: the number base pairs to extend to upstream and downstream. Default 5.\n");
	  fprintf(stderr,"       --noGapInRead: when use every nucleotide in the extended read instead of just the nucleotides on two ends\n");
	  fprintf(stderr,"       --gc_bin: bin the data and report gc and read count in the bins; Notice that if the option -e is specified, this option will be ignored\n");
	  fprintf(stderr,"       --map_bin: if specified, bin the data as bins with equal number of uniquely mappable genomic locations\n");
	  fprintf(stderr,"                  only valid if --gc_bin is specified and this option assume the input file has all unique mappable positions of the chromosome under consideration\n");
	  fprintf(stderr,"       --NoNegBinom: Do not use negative binomial to estimate the variance;\n");
	  fprintf(stderr,"                     If this is not specified, the last row of -e <string> will be used as estimate of theta parameter in the Negative Binomial model\n");
	  fprintf(stderr,"       --NoHeader: Do not print header to the output\n");
	  fprintf(stderr,"       --refine <string>: the file that stores the parameters in the refine step\n");
          return;
        }



static DPargs  option_assign(int argc, char **argv)
        { DPargs args;
          int c,option_index;
	  char *poisfile=NULL,*arg_refine=NULL;
	  FILE *refineFile = NULL;
	  double *x;
	  int nrow,ncol;
          static struct option long_options[] =
                        { {"chrom",required_argument,0,0}, 
			  {"gc_bin",no_argument,0,1},
			  {"map_bin",no_argument,0,2},
			  {"NoNegBinom",no_argument,0,4},
			  {"refine",required_argument,0,6},
			  {"NoHeader",no_argument,0,7},
			  {"uds", required_argument,0,8},
			  {"noGapInRead",no_argument,0,9},
                          {0,0,0,0}
                        };


          args.fa_file_name = NULL;
          args.readcnt_file_name = NULL;
          args.output_name = NULL;
          args.chrom_name=NULL;
          args.fragment_size = 300;
	  args.readLen = 50;
	  args.extd_bp = 5;
	  args.input_pois = NULL;
	  args.estimate = 0;
	  args.bin_size = 100;
	  args.gc_bin = 0;
	  args.negBinom = 1;
	  args.map_bin = 0;
	  args.fa_file = NULL;
	  args.rdcnt_file = NULL;
	  args.output = NULL;
	  args.inverse_theta = 0;
	  args.nogapInRead = 0;

	  args.refine = 0;
	  args.rf_gc_coef = NULL;
	  args.size_rf_gc_coef = 0;

          while((c=getopt_long(argc,argv,"g:i:s:e:b:l:o:h",long_options,&option_index))!=-1){
                switch(c){
                        case 0:
                                args.chrom_name = strdup(optarg);
				break;
			case 1:
				args.gc_bin = 1;
                                break;
			case 2:
				args.map_bin = 1;
				break;
			case 4: 
				args.negBinom = 0;
				break;
			case 6:
				arg_refine = strdup(optarg);
				args.refine = 1;
				break;
			case 7:
				no_header=1;
				break;
			case 8:
				args.extd_bp = atoi(optarg);
				if(args.extd_bp<0){fprintf(stderr,"The option --uds only takes nonnegative values\n"); exit(1);}
				break;
			case 9:
				args.nogapInRead = 1;
				break;
			case 'i':
				args.readcnt_file_name = strdup(optarg);
				break;
                        case 's':
                                args.fragment_size = atoi(optarg);
				if(args.fragment_size<=0) {fprintf(stderr,"Error fragment size (=%d) must be positive\n",args.fragment_size);exit(1);}
                                break;
                        case 'l':
                                args.readLen = atoi(optarg);
                                if(args.readLen<0) {fprintf(stderr,"The read Length must be nonnegative.\n");exit(1);}
                                break;
                        case 'o':
                                args.output_name = strdup(optarg);
                                break;
			case 'e':
				poisfile = strdup(optarg);
				args.estimate = 1;
				break;
			case 'b':
				args.bin_size = atoi(optarg);
				if(args.bin_size<=0) {fprintf(stderr,"The bin size must be postive\n");exit(1);}
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

        args.fa_file_name = strdup(argv[optind]);

	if(poisfile!=NULL){
		args.input_pois = fopen(poisfile,"r");
		if(args.input_pois==NULL){
			fprintf(stderr,"Failed to open the file: %s\n",poisfile);
			exit(1);
			}
		free(poisfile); poisfile=NULL;
		}

        if(args.readcnt_file_name!=NULL){
                args.rdcnt_file = fopen(args.readcnt_file_name,"r");
                if(args.rdcnt_file==NULL){
                        fprintf(stderr,"fopen %s: %s\n", args.readcnt_file_name,strerror(errno));
                        exit(1);
                        }
                }else{
                args.rdcnt_file = stdin;
                }


        args.fa_file = fopen(args.fa_file_name,"r");
        if(args.fa_file==NULL)
                { fprintf(stderr,"fopen %s: %s\n", args.fa_file_name,strerror(errno));
                  exit(1);
                }

        if(args.output_name!=NULL){
                args.output = fopen(args.output_name,"w");
                if(args.output==NULL){
                        fprintf(stderr,"fopen %s: %s\n",args.output_name,strerror(errno));
                        exit(1);
                        }
                }
        else args.output = stdout;


	if(args.refine==1){
		int i;
		refineFile = fopen(arg_refine,"r");
		if(refineFile==NULL){
			fprintf(stderr,"Failed to open the file: %s\n",arg_refine);
			exit(1);
			}
		x = read_table(refineFile,&nrow,&ncol,-1,0);
		if(nrow!=2*args.fragment_size+2){
			fprintf(stderr,"Expected to get %d rows numeric values in %s, but get %d rows\n",2*args.fragment_size+2,arg_refine,nrow);
			exit(1);
			}

		args.size_rf_gc_coef = nrow;
		args.rf_gc_coef = (double *)malloc(sizeof(double)*(nrow+2));
		for(i=0;i<nrow;i++){
			args.rf_gc_coef[i] = x[i*ncol];
			}

		}

        if(args.readLen<0||args.readLen>args.fragment_size){
                fprintf(stderr,"<Read length> must be positive and smaller than the <Fragment size>\n");
                }

	if(args.nogapInRead==0&& 2*args.extd_bp > args.readLen){
		fprintf(stderr,"The value of option \'uds\' must be no larger than a half of readLenth if \'--noGapInRead is specified\'");
		}
	if(args.nogapInRead==0){
		args.NumNuclotide = 4*args.extd_bp;
		}else{
		args.NumNuclotide = 2*args.extd_bp + args.readLen;
		}

        fprintf(stderr,"Fragment length = %d\n",args.fragment_size);
        fprintf(stderr,"Read Length = %d\n",args.readLen);


        return args;

        }




char *read_fa(FILE *fa_file,int *len,char **chrom)
	{ char *dna,chr,*ll;
	  int min_size = 260000000;
	  int i=0;

	  dna = (char *)malloc(sizeof(char)*min_size);
	  if(dna==NULL) { fprintf(stderr,"Error: memory allocation failed\n");exit(1);}

	  chr=getc(fa_file);
	  if(chr=='>')
	  	{ ll = lr_read(fa_file); /*this line is the chromosome information*/
		  *chrom = strdup(ll);
		}
	  else { ungetc(chr,fa_file);
		 *chrom =(char*)malloc(sizeof(char)*10);
		 if(*chrom ==NULL) {fprintf(stderr,"Error: memory allocation failed\n");exit(1);}
		 *chrom[0] = '\0';
		}

	  i = 0;
	  while((chr=getc(fa_file))!=EOF&&chr!='>')
		{ if(chr!='\n'){
			dna[i] = chr;
			i++;
			}
		  if(i+10>=min_size){
			min_size += 10000000;
			dna = (char *) realloc(dna,sizeof(char)*min_size);
			if(dna==NULL)  { fprintf(stderr,"Error: memory allocation failed\n");exit(1);}
			}
		}
	  if(chr=='>') ungetc(chr,fa_file);

	  *len = i;
	
	  return dna;
	}



void print_header(DPargs args)
	{ int i;
	  FILE *output = args.output;
	  int size = args.readLen + 2*args.extd_bp;

	  if(output==NULL) return;
	  fprintf(output,"pos\tcount");
	  fprintf(output,"\tGC");

	  if(args.nogapInRead==1){	
		  for(i=0;i<size;i++){
			fprintf(output,"\tG%d\tC%d\tT%d",i+1,i+1,i+1);
			}
		}else{
		  for(i=0;i<2*args.extd_bp;i++){
			fprintf(output,"\tG_Start%d\tC_Start%d\tT_Start%d",i+1,i+1,i+1); /*[-args.extd_bp, args.extd_bp-1] nucleotide*/
			}
		  for(i=0;i<2*args.extd_bp;i++){
			fprintf(output,"\tG_End%d\tC_End%d\tT_End%d",i+1,i+1,i+1); /*[readLen-args.extd_bp, readLen+args.extd_bp-1] nucleotide*/
			}
		}
	  fprintf(output,"\n");
	  return;
	}


int count_nuc(char *dna, int length, int start_in, int end_in, char c)
	{ int cnt = 0,i;
	  int start, end;
	  int ll = length-1;
	  start = start_in >0 ? start_in : 0; // modified here in Version 0.1.3
	  end = end_in < ll ? end_in: ll;
	  for(i=start;i<=end;i++){
		if(toupper(dna[i])==toupper(c)) cnt++;
		}
	  return cnt;
	}



void print_oneline(DPargs args, char *dna, int dna_len, int pos_in,int read_count)
	{int gc, g, c, t,i,pos, k;
	 double gc_fraction,tmp;
	 FILE *output = args.output;
	 int fragment_size=args.fragment_size, read_len = args.readLen;
	 int size = args.readLen + 2*args.extd_bp;

	pos = pos_in;

	if(fragment_size<=0||read_len<0){return;}
	if(pos<read_len||pos+read_len>=dna_len||pos<fragment_size||pos+fragment_size>=dna_len){return;}

	gc = count_nuc(dna,dna_len,pos-fragment_size,pos+fragment_size,'G');
	gc += count_nuc(dna,dna_len,pos-fragment_size,pos+fragment_size,'C');
	
	gc_fraction = ((double) gc)/((double) 2*fragment_size+1);

	fprintf(output,"%d\t%d",pos_in,read_count);

	tmp = gc_fraction;
	fprintf(output,"\t%g",tmp);

	if(args.nogapInRead==1){ /*do not put gap in the read for model training and prediction*/
		for(i=0;i<size;i++){
			k = i - args.extd_bp;
			g = count_nuc(dna,dna_len,pos+k,pos+k,'g');
			c = count_nuc(dna,dna_len,pos+k,pos+k,'c');
			t = count_nuc(dna,dna_len,pos+k,pos+k,'t');
			fprintf(output,"\t%d\t%d\t%d",g,c,t);
			}
		}else{
		for(i=0;i<2*args.extd_bp;i++){ /*first end*/
			k = i - args.extd_bp;
                        g = count_nuc(dna,dna_len,pos+k,pos+k,'g');
                        c = count_nuc(dna,dna_len,pos+k,pos+k,'c');
                        t = count_nuc(dna,dna_len,pos+k,pos+k,'t');
                        fprintf(output,"\t%d\t%d\t%d",g,c,t);
			}

		for(i=0;i<2*args.extd_bp;i++){ /*second end*/
                        k = i + args.readLen - args.extd_bp;
                        g = count_nuc(dna,dna_len,pos+k,pos+k,'g');
                        c = count_nuc(dna,dna_len,pos+k,pos+k,'c');
                        t = count_nuc(dna,dna_len,pos+k,pos+k,'t');
                        fprintf(output,"\t%d\t%d\t%d",g,c,t);
			}
		}

	fprintf(output,"\n");
	return;
	}

static int update_gc(char *dna, int dna_len, int pos, int pos_pre, int fragment_size, int gc_pre)
{	int  gc;

	if(pos-fragment_size > pos_pre+fragment_size || pos_pre < 0){
		gc = count_nuc(dna,dna_len,pos-fragment_size,pos+fragment_size,'G');
		gc += count_nuc(dna,dna_len,pos-fragment_size,pos+fragment_size,'C');
		}
	else {
		gc = gc_pre;
		gc -= count_nuc(dna,dna_len,pos_pre-fragment_size,pos-fragment_size-1,'G');
		gc -= count_nuc(dna,dna_len,pos_pre-fragment_size,pos-fragment_size-1,'C');
		gc += count_nuc(dna,dna_len,pos_pre+fragment_size+1,pos+fragment_size,'G');
		gc += count_nuc(dna,dna_len,pos_pre+fragment_size+1,pos+fragment_size,'C');
		}

	return gc;
}


double calculate_bias(DPargs args, char *dna, int dna_len, int pos_in, double *param_est, int nrow_param, int ncol_param, double *gc_bias)
{	double bias = 0.0;
	int gc, g, c, t, i,pos, k;
        double gc_fraction;
	int fragment_size = args.fragment_size;
	int size = args.readLen + 2*args.extd_bp;
	static int pos_pre=-1, gc_pre = 0;
	//int gc_tmp;

	pos = pos_in;
	//if(pos<read_len||pos+read_len>=dna_len||pos<fragment_size||pos+fragment_size>=dna_len){return 0.0;}

	gc = update_gc(dna, dna_len,pos,pos_pre,fragment_size,gc_pre);
	if(gc >= 2*fragment_size+2){fprintf(stderr,"Error, gc = %d should be less than 2*fragmentSize+2 = %d", gc, 2*fragment_size+2); exit(1);} /*Just in case, later I should delete this*/

	//gc_tmp = count_nuc(dna,dna_len,pos-fragment_size,pos+fragment_size,'G');
	//gc_tmp += count_nuc(dna,dna_len,pos-fragment_size,pos+fragment_size,'C');

	//if(gc!=gc_tmp) {fprintf(stderr,"Error, updated gc = %d, newly calculated gc = %d, pos = %d, pos_pre = %d\n", gc, gc_tmp, pos, pos_pre); exit(1);}

	pos_pre = pos; gc_pre = gc;

	gc_fraction = ((double) gc)/((double) 2*fragment_size+1);
	
	bias = param_est[ncol_param*(gc+3*args.NumNuclotide)];
	*gc_bias = exp(bias);	

	if(args.nogapInRead==1){ /*do not put gap in the read for model training and prediction*/
		for(i=0;i<size;i++){
			k = i - args.extd_bp;
       	        	g = count_nuc(dna,dna_len,pos+k,pos+k,'g');
                	c = count_nuc(dna,dna_len,pos+k,pos+k,'c');
                	t = count_nuc(dna,dna_len,pos+k,pos+k,'t');
		
			bias += g*param_est[(3*i)*ncol_param];
			bias += c*param_est[(3*i+1)*ncol_param];
			bias += t*param_est[(3*i+2)*ncol_param];
			}
		}else{
		int j=0;/*j is to index the parameter*/
		for(i=0;i<2*args.extd_bp;i++){ /*First end*/
                        k = i - args.extd_bp;
                        g = count_nuc(dna,dna_len,pos+k,pos+k,'g');
                        c = count_nuc(dna,dna_len,pos+k,pos+k,'c');
                        t = count_nuc(dna,dna_len,pos+k,pos+k,'t');

                        bias += g*param_est[(3*j)*ncol_param];
                        bias += c*param_est[(3*j+1)*ncol_param];
                        bias += t*param_est[(3*j+2)*ncol_param];

			j++;
                        }

                for(i=0;i<2*args.extd_bp;i++){ /*Second end*/
                        k = i + args.readLen - args.extd_bp;
                        g = count_nuc(dna,dna_len,pos+k,pos+k,'g');
                        c = count_nuc(dna,dna_len,pos+k,pos+k,'c');
                        t = count_nuc(dna,dna_len,pos+k,pos+k,'t');

                        bias += g*param_est[(3*j)*ncol_param];
                        bias += c*param_est[(3*j+1)*ncol_param];
                        bias += t*param_est[(3*j+2)*ncol_param];

                        j++;
                        }
		}

	if(args.refine==1){
		bias += args.rf_gc_coef[gc];;
		}

	return exp(bias);
	}



void fprint_gc_bin(FILE *output, char*dna, int dna_len, double *rdcnt, int nrow_rdcnt,  int bin_size, int map_bin, char *chrom_name_in)
{	int i,j,end,start, pos,pos_pre;
	double gc_cnt,n_cnt, rdcnt_tmp;
	//double dist1, dist2;
	char *chrom_name=NULL;

	if(chrom_name_in!=NULL) {
		int len = strlen(chrom_name_in);
		chrom_name = (char *)malloc(sizeof(char)*(len+10));
		chrom_name = strcpy(chrom_name,chrom_name_in);
		chrom_name[len] = '\t';
		chrom_name[len+1] = '\0';
		}
	else{ 
		chrom_name = (char *) malloc(sizeof(char)*10);
		chrom_name[0] = '\0';
		}

	if(no_header!=1) {
		if(chrom_name_in!=NULL){fprintf(output,"chrom\t");}
		fprintf(output,"start\tend\tread_count\tgc");
		fprintf(output,"\n");
		}
	if(map_bin==1){
		int num_mapbin = nrow_rdcnt/bin_size;
		//for(i=0;i<nrow_rdcnt/bin_size+1;i++){
		for(i=0;i<num_mapbin;i++){
			start = i*bin_size;
			end = (i+1)*bin_size-1;
			//end = end<nrow_rdcnt-1 ? end : nrow_rdcnt-1;
			if(i==num_mapbin-1) end = nrow_rdcnt-1;
			pos_pre = (int) rdcnt[start*2];
			pos = (int) rdcnt[end*2];
			gc_cnt = (double) count_nuc(dna, dna_len,pos_pre-1,pos-1,'g');
			gc_cnt += (double) count_nuc(dna, dna_len,pos_pre-1,pos-1,'c');
			n_cnt = (double) count_nuc(dna, dna_len,pos_pre-1,pos-1,'n');;
			rdcnt_tmp = 0.0;
			for(j=start;j<=end;j++){
				rdcnt_tmp +=  rdcnt[j*2+1];
				}
			if(end>=start&&n_cnt<pos-pos_pre+1) {
				fprintf(output,"%s%d\t%d\t%.0f\t%g",chrom_name,pos_pre,pos,rdcnt_tmp,gc_cnt/(pos-pos_pre+1-n_cnt));
				fprintf(output,"\n");
				}
			}
		}else{
		gc_cnt = 0.0;
		n_cnt = 0.0;
		rdcnt_tmp = 0.0;
		start = (int) floor((rdcnt[0]-1)/bin_size)*bin_size +1;
		end = start + bin_size -1;
		for(i=0;i<nrow_rdcnt;i++){
			pos = (int) rdcnt[2*i];
			rdcnt_tmp += rdcnt[2*i+1];
			if(pos > end||i==nrow_rdcnt-1){
				gc_cnt = (double) count_nuc(dna, dna_len,start-1,end-1,'g');
				gc_cnt += (double) count_nuc(dna, dna_len,start-1,end-1,'c');
				n_cnt = (double) count_nuc(dna, dna_len,start-1,end-1,'n');
				if((int) n_cnt < bin_size ){
					fprintf(output,"%s%d\t%d\t%.0f\t%g",chrom_name,start, end, rdcnt_tmp, gc_cnt/(bin_size-n_cnt));
					fprintf(output,"\n");
					}
				rdcnt_tmp = 0.0;
				start = bin_size*(pos/bin_size)+1;
				end = start + bin_size-1;
				}
			}
		}
	return ;
}


void fprintf_predicted(DPargs args, double *rdcnt, int nrow_rdcnt,  char *dna, int dna_len, double *param_est, int nrow_param, int ncol_param)
{	int bin_start, bin_end, count_sum,pos,pos_pre,i,read_count,j;
	double  cnt_exp, gc_bias_tmp, gc_cnt , n_cnt,var,mean_tmp;


        count_sum = 0;
        cnt_exp = 0.0;
	var = 0.0;

	if(no_header!=1){
		fprintf(args.output,"start\tend\tobs\texpected\tvar");
		if(args.gc_bin==1){fprintf(args.output,"\tgc");}
		fprintf(args.output,"\n");
		}

	if(args.map_bin!=1){
		if(nrow_rdcnt>0){
			bin_start = (int) floor((rdcnt[0]-1)/args.bin_size)*args.bin_size +1; 
			bin_end = bin_start+args.bin_size-1;
			}
		for(i=0;i<nrow_rdcnt;i++){
			pos = (int) rdcnt[2*i];
			read_count = (int) rdcnt[2*i+1];

			if(pos<0||pos>dna_len){
				fprintf(stderr,"No such position (%d) in fa file %s\n",pos,args.fa_file_name);
				exit(1);
				}
			if(read_count<0){
				fprintf(stderr,"Error, negative read count (%d), line %d in %s\n",read_count,i,args.readcnt_file_name);
				exit(1);
				}

			mean_tmp =  calculate_bias(args, dna, dna_len, pos, param_est, nrow_param, ncol_param, &gc_bias_tmp);

			if(pos<=bin_end && i!=nrow_rdcnt-1){
				cnt_exp += mean_tmp;
				var += mean_tmp + args.inverse_theta*mean_tmp*mean_tmp;
				count_sum += read_count;
				}else{
				if(args.gc_bin==1){
					gc_cnt = (double)  count_nuc(dna, dna_len,bin_start-1,bin_end-1,'g');
					gc_cnt += (double) count_nuc(dna, dna_len,bin_start-1,bin_end-1,'c');
					n_cnt = (double)  count_nuc(dna, dna_len,bin_start-1,bin_end-1,'n');
					if(n_cnt < args.bin_size) {
						fprintf(args.output,"%d\t%d\t%d\t%g\t%g",bin_start,bin_end,count_sum,cnt_exp,var);
						fprintf(args.output,"\t%g\n",gc_cnt/(args.bin_size-n_cnt));
						}
					}else{
					fprintf(args.output,"%d\t%d\t%d\t%g\t%g\n",bin_start,bin_end,count_sum,cnt_exp,var);
					}
				cnt_exp = mean_tmp;
				var = mean_tmp + args.inverse_theta*mean_tmp*mean_tmp;
				count_sum = read_count;
				bin_start = (pos/args.bin_size)*args.bin_size+1;
				bin_end = bin_start+args.bin_size-1;
				}
			}
	}else{
        int num_mapbin = nrow_rdcnt/args.bin_size;
        //for(i=0;i<nrow_rdcnt/args.bin_size+1;i++){
        for(i=0;i<num_mapbin;i++){
                bin_start = i*args.bin_size;
                bin_end = (i+1)*args.bin_size-1;
                //bin_end = (bin_end < nrow_rdcnt-1 ? bin_end : nrow_rdcnt-1);
                if(i==num_mapbin-1) bin_end = nrow_rdcnt-1;
                pos_pre = (int) rdcnt[bin_start*2];
                pos = (int) rdcnt[bin_end*2];
		count_sum = 0;
		cnt_exp = 0.0;
		//gc_bias = 0.0;
		//cnt_normal = 0.0;
		var = 0.0;
		for(j=bin_start;j<=bin_end;j++){
			mean_tmp =  calculate_bias(args, dna, dna_len, (int) rdcnt[2*j], param_est, nrow_param, ncol_param, &gc_bias_tmp);
			//gc_bias += gc_bias_tmp;
			cnt_exp += mean_tmp;
			var += mean_tmp + args.inverse_theta*mean_tmp*mean_tmp;
			count_sum += (int) rdcnt[j*2+1];
			}
		if(bin_end>=bin_start){
			//cnt_normal = ((double) count_sum)/(cnt_exp);
			if(args.gc_bin==1){
                		gc_cnt = (double) count_nuc(dna, dna_len,pos_pre-1,pos-1,'g');
		                gc_cnt += (double) count_nuc(dna, dna_len,pos_pre-1,pos-1,'c');
		                n_cnt = (double) count_nuc(dna, dna_len,pos_pre-1,pos-1,'n');
				if(n_cnt<pos-pos_pre+1) {
					fprintf(args.output,"%d\t%d\t%d\t%g\t%g",pos_pre,pos,count_sum,cnt_exp,var);
					fprintf(args.output,"\t%g\n",gc_cnt/(pos-pos_pre+1-n_cnt));
					}
				}else{
				fprintf(args.output,"%d\t%d\t%d\t%g\t%g\n",pos_pre,pos,count_sum,cnt_exp,var);
				}
			}

		}
	}
}


