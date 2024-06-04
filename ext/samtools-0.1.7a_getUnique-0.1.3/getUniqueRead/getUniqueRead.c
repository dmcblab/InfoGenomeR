
#include "getUniqueRead.h"
#include "../bam.h"

#define Bowtie 1
#define BWA 2

typedef struct ReadPosOutput_t{
	char *chrom;
	char *outputname;
	FILE *output;
	int num_reads; /*number of reads reported so far for this chromosome*/
	} ReadPosOutput;

typedef struct ReadPosOutput_allchrom_t{
	ReadPosOutput * cout;
	int size;
	} RdOutAC;

static char *outputPrefix = NULL;
static int ALIGNER = BWA;
static int ChromNameReport = 0;
static int StrandReport = 0;
static int StrandSep = 0; /*Separate negative reads and positive reads to separate files? (0 no, 1 yes)*/
static int minReadLen = 0;
static int maxReadLen = 100000000;

static int *readLength=NULL;
static int maxReadLen_S = 0;


void summarizeReadLength(bam1_t *b)
{	int k = 0, i=0, mRL;
	k = b->core.l_qseq;
	//printf("k=%d\n",k);
	if(readLength==NULL || k +10> maxReadLen_S){
		mRL = maxReadLen_S;
		maxReadLen_S = k+10;
		readLength = (int *) realloc(readLength,sizeof(int)*(maxReadLen_S+10));
		for(i=mRL;i<maxReadLen_S;i++){
			readLength[i]=0;
			}
		}
	readLength[k-1]++;

  return;
}

void reportReadLength(FILE *out)
{	int i;
	for(i=0;i<maxReadLen_S;i++){
		if(i==0) fprintf(out,"readLength\tcount\n");
		if(readLength[i]>0){
			fprintf(out,"%d\t%d\n",i+1,readLength[i]);
			}
		}

}

//static int cmp_ReadPosOutput(void *a1, void *a2);
static int is_unique_mapped(bam1_t *b, int aligner); /*In this version, this function will just test if a read is mapped*/
static RdOutAC *create_Readoutput(bam_header_t *in, char *prefix);
static int destroy_Readoutput(RdOutAC * out); /*return the total number of reads reported so far*/
static int report_one_read(RdOutAC *out, bam1_t *b);

static RdOutAC *create_Readoutput(bam_header_t *in, char *prefix)
{	RdOutAC *ot=malloc(sizeof(RdOutAC));
	int i;
	ot->size = (in->n_targets)*(StrandSep+1);
	ot->cout = malloc(sizeof(ReadPosOutput)*(ot->size+10));
	if(ot->cout==NULL) {fprintf(stderr,"Error in create_Readoutput: memory allocation failed\n");exit(1);}
	
	if(StrandSep==0){
		for(i=0;i<ot->size;i++){
			ot->cout[i].chrom = strdup(in->target_name[i]);
			ot->cout[i].outputname = malloc(sizeof(char)*(strlen(prefix)+strlen(ot->cout[i].chrom)+10));
			if(ot->cout[i].outputname==NULL) {fprintf(stderr,"Error in create_Readoutput: memory allocation failed\n");exit(1);}
			ot->cout[i].outputname[0] = '\0';
			strcat(ot->cout[i].outputname,prefix);
			strcat(ot->cout[i].outputname,ot->cout[i].chrom);
			strcat(ot->cout[i].outputname,".seq");

			ot->cout[i].output=NULL; 
			ot->cout[i].output = fopen(ot->cout[i].outputname,"w");
			if(ot->cout[i].output==NULL){
				fprintf(stderr,"Error, failed to create/open the file %s\n",ot->cout[i].outputname);
				exit(1);
				}

			ot->cout[i].num_reads=0;
			}
		}else{
                for(i=0;i<ot->size/2;i++){
                        ot->cout[2*i].chrom = strdup(in->target_name[i]);
                        ot->cout[2*i].outputname = malloc(sizeof(char)*(strlen(prefix)+strlen(ot->cout[2*i].chrom)+20));
                        if(ot->cout[2*i].outputname==NULL) {fprintf(stderr,"Error in create_Readoutput: memory allocation failed\n");exit(1);}
                        ot->cout[2*i].outputname[0] = '\0';
                        strcat(ot->cout[2*i].outputname,prefix);
                        strcat(ot->cout[2*i].outputname,ot->cout[2*i].chrom);
			strcat(ot->cout[2*i].outputname,"_pos");
                        strcat(ot->cout[2*i].outputname,".seq");

                        ot->cout[2*i].output = fopen(ot->cout[2*i].outputname,"w");
                        if(ot->cout[2*i].output==NULL){
                                fprintf(stderr,"Error, failed to create/open the file %s\n",ot->cout[2*i].outputname);
                                exit(1);
                                }

                        ot->cout[2*i].num_reads=0;

                        ot->cout[2*i+1].chrom = strdup(in->target_name[i]);
                        ot->cout[2*i+1].outputname = malloc(sizeof(char)*(strlen(prefix)+strlen(ot->cout[2*i].chrom)+20));
                        if(ot->cout[2*i+1].outputname==NULL) {fprintf(stderr,"Error in create_Readoutput: memory allocation failed\n");exit(1);}
                        ot->cout[2*i+1].outputname[0] = '\0';
                        strcat(ot->cout[2*i+1].outputname,prefix);
                        strcat(ot->cout[2*i+1].outputname,ot->cout[2*i].chrom);
                        strcat(ot->cout[2*i+1].outputname,"_neg");
                        strcat(ot->cout[2*i+1].outputname,".seq");

                        ot->cout[2*i+1].output = fopen(ot->cout[2*i+1].outputname,"w");
                        if(ot->cout[2*i+1].output==NULL){
                                fprintf(stderr,"Error, failed to create/open the file %s\n",ot->cout[2*i+1].outputname);
                                exit(1);
                                }

                        ot->cout[2*i+1].num_reads=0;

                        }


		}

	return ot;
}

static int destroy_Readoutput(RdOutAC *out)
{	int i, n=0;
	if(out==NULL) return 0;

	for(i=0;i<out->size;i++){
		n += out->cout[i].num_reads;
		//fprintf(stderr,"chromsome %s: %d reads\n",out->cout[i].chrom,out->cout[i].num_reads);
		free(out->cout[i].chrom);
		free(out->cout[i].outputname);
		fclose(out->cout[i].output);
		}
	free(out->cout);
	out->cout=NULL;
	out->size = 0;
	
	free(out);
	return n;
}


int get_Unique_args(char *args_in)
{	char *args = strdup(args_in), *substr=NULL;
 	if(args==NULL){
		fprintf(stderr,"argument is required for option U\n");
		exit(1);
		}

	substr = strtok(args,",");
	if(strcmp(substr,"BWA")==0){
		ALIGNER = BWA;
		//fprintf(stderr,"Aliger is BWA\n");
		}else if(strcmp(substr,"Bowtie")==0){
		ALIGNER = Bowtie;
		//fprintf(stderr,"Aliger is Bowtie\n");
		}else{
		fprintf(stderr,"Aligner can only be BWA or Bowtie\n");
		exit(1);
		}
	substr = strtok(NULL,",");
	if(substr==NULL){
		fprintf(stderr, "Option U should have argument of the format\n<Aligner,OutputPrefix,ChomNameReport?,StrandReport?> or <Aligner,OutputPrefix,ChomNameReport?,StrandReport?,minLen,MaxLen>\n");
		exit(1);
		}
	outputPrefix = strdup(substr);

	substr = strtok(NULL,",");
        if(substr==NULL){
                fprintf(stderr, "Option U should have argument of the format\n<Aligner,OutputPrefix,ChomNameReport?,StrandReport?> or <Aligner,OutputPrefix,ChomNameReport?,StrandReport?,minLen,MaxLen>\n");
                exit(1);
                }
	if(strcmp(substr,"N")!=0&&strcmp(substr,"Y")!=0){
		fprintf(stderr,"The 3rd argumennt for the option U must be 'Y' or 'N'");
		exit(1);
		}else if(strcmp(substr,"N")==0){ChromNameReport=0;}
		else if(strcmp(substr,"Y")==0){ChromNameReport=1;}
	
	substr = strtok(NULL,",");
        if(substr==NULL){
                fprintf(stderr, "Option U should have argument of the format\n<Aligner,OutputPrefix,ChomNameReport?,StrandReport?> or <Aligner,OutputPrefix,ChomNameReport?,StrandReport?,minLen,MaxLen>\n");
                exit(1);
                }
        if(strcmp(substr,"N")!=0&&strcmp(substr,"Y")!=0&&strcmp(substr,"S")!=0){
                fprintf(stderr,"The 4th argumennt for the option U must be 'Y' or 'N', or 'S'");
		exit(1);
                }else if(strcmp(substr,"N")==0){StrandReport=0;}
                else if(strcmp(substr,"Y")==0){StrandReport=1;}
		else if(strcmp(substr,"S")==0){StrandSep=1;}

	substr = strtok(NULL,",");
	if(substr!=NULL){
		minReadLen = atoi(substr);
		if(minReadLen<0) fprintf(stderr,"minLen=%d is negative\n",minReadLen);

		substr = strtok(NULL,",");
	        if(substr==NULL){
        	        fprintf(stderr, "Option U should have argument of the format\n<Aligner,OutputPrefix,ChomNameReport?,StrandReport?> or <Aligner,OutputPrefix,ChomNameReport?,StrandReport?,minLen,MaxLen>\n");
               		exit(1);
                	}
		maxReadLen = atoi(substr);
		if(maxReadLen<minReadLen) fprintf(stderr,"maxLen=%d is less than minLen=%d\n",maxReadLen, minReadLen);
		}



	return 0;
}
/*
static int cmp_ReadPosOutput(void *a1, void *a2)
{	ReadPosOutput tmp1, tmp2;
	tmp1 = *(ReadPosOutput *) a1;
	tmp2 = *(ReadPosOutput *) a2;

	if(strcmp(tmp1.chrom,tmp2.chrom)<0) return -1;
	else if (strcmp(tmp1.chrom,tmp2.chrom)>0) return 1;
	else return 0;
}*/


static int is_unique_mapped(bam1_t *b, int aligner)
{	int is_unmapped;
	//int num_hit;
	//char XTtag;
	//uint8_t *tag_p;


	if(b->core.l_qseq < minReadLen || b->core.l_qseq > maxReadLen) return 0; /*Does not satisfy the read length constraint, return false (0)*/

	is_unmapped = ((b->core.flag) & 0x0004)>>2;

	if(is_unmapped==0){
		//if(aligner==Bowtie)
		//	{ tag_p = bam_aux_get(b,"XM");
		//	  num_hit = bam_aux2i(tag_p);
		//	  if(num_hit<=1) return 1; else return 0;
		//	}
		//else if(aligner==BWA){
		//	tag_p = bam_aux_get(b,"XT");
		//	XTtag = bam_aux2A(tag_p);
		//	tag_p = bam_aux_get(b,"X1");
		//	num_hit = bam_aux2i(tag_p);

			//if(XTtag!='U'||num_hit>1) {fprintf(stderr,"name=%s\tXTtag=%c\tnum_hit=%d\n",bam1_qname(b),XTtag,num_hit);getchar();}

		//	if(XTtag=='U'&&num_hit<=0) return 1; 
		//	else return 0;

		//	}
		return 1;
	}

	return 0;
}


static int report_one_read(RdOutAC *out, bam1_t *b)
{	int unique;
	int strand;
	FILE *fout=NULL;
	int index;
	unique = is_unique_mapped(b,ALIGNER);
	if(unique==0) return 0; /*b is not unqiuely mapped, do not report*/

	strand = ((b->core.flag) & 0x0010)>>4;
	if(StrandSep==1){
		fout = out->cout[b->core.tid*2+strand].output;
		index = b->core.tid*2+strand;
		}else{
		fout = out->cout[b->core.tid].output;
		index = b->core.tid;
		}
	if(fout==NULL) {fprintf(stderr,"ERROR: fout is NULL\n");}

	if(ChromNameReport==1) fprintf(fout,"%s\t",out->cout[index].chrom);
	fprintf(fout,"%d",b->core.pos);
	if(StrandReport==1) {
		fprintf(fout,"\t%d",strand);
		}
	fprintf(fout,"\n");

	out->cout[index].num_reads++;
	return 1;
}

void report_unique_read(bam_header_t *in,bam1_t *b)
{	static RdOutAC *out= NULL;
	int num_unique=0;
	
	if(out==NULL) out = create_Readoutput(in, outputPrefix);

	if(b!=NULL) report_one_read(out,b);


	if(in==NULL&&b==NULL&&out!=NULL) {
		num_unique=destroy_Readoutput(out);
		fprintf(stderr,"Reported %d unqiue mapped reads\n",num_unique);
		}
	return ;
}
