/* Written by Ruibin Xi, April 16 2010
 *
 * */
#ifndef SINGULARITY_RM_XI
#define SINGULARITY_RM_XI

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct COUNT_LST_NODE_t {
		int num_reads; /*number of reads at a position*/
		int num_pos; /*how many position in this window that has num_reads number of reads*/
		struct COUNT_LST_NODE_t *next, *pre;
	} *COUNT_LST_NODE; 
/*a node of list that represents the distribution of read count at the positions in the window*/

typedef struct COUNT_LST_t{
	 COUNT_LST_NODE head;
	 COUNT_LST_NODE tail;
	 int w; /*the size of the window*/
	} *COUNT_LST;

/*the list COUNT_LIST will be ordered decreasingly by its element num_reads (from head to tail)*/

typedef struct POS_LST_NODE_t{
		int pos;
		int read_cnt;
		struct POS_LST_NODE_t *next, *pre;
		COUNT_LST_NODE cnt_node;
	} *POS_LST_NODE;

typedef struct POS_LST_t{
		POS_LST_NODE head;
		POS_LST_NODE tail;
		int len; /*the length of this list, should be equal to window size*/
	} *POS_LST;


typedef struct LST_t{
		POS_LST pos_lst;
		COUNT_LST cnt_lst;
	} *LST;

typedef struct Quantile_t{
		double p; /*the pth  quantile*/
		int w; /*the window size*/
		COUNT_LST_NODE cnt_node; /*the count node that corresponding to the pth quantile*/
		int k;/*number of genomic positions after this node(exclusive)*/		
	} *Quantile;

/***************************************************************************
 ***************************************************************************
 ************** functions for count list ***********************************
 ***************************************************************************/

COUNT_LST COUNT_LST_create(int win_size); /*create a COUNT_LST*/

COUNT_LST_NODE  insert_CNT(COUNT_LST ll, int read_count);
/* Given a read count, 
 * 1. if there is a node in the list whose element num_reads equals to 
 *    read_count, then increase num_pos of this node by 1.
 * 2. if there is no such a node, create a new node with its num_reads 
 *    as read_count and its num_pos as 1 while keeping the order of the list
 * 3. if ll==NULL, report an error;
 *
 * The returned value is the node that was UPDATED, but not the head of the list
 * */

void delete_CNT(COUNT_LST ll, COUNT_LST_NODE node);
/* similar to insert_RDL, but this function is for deletion.
 * */



/***************************************************************************/
/***************************************************************************/
/*************** functions for position list *******************************/
/***************************************************************************/


POS_LST POS_LST_create();

void push_POS(POS_LST ll, int pos, int read_cnt); /*push a node on the top of the list*/
void del_tail_POS(POS_LST ll); /*delete the tail of the list*/



/********************************************
 ********************************************
 ****** functions for quantile calculation***
 ********************************************/

LST LST_create(int w); /*w is the window size*/



void lst_initialize(int *read_pos, int num_pos, LST lst, int win_size);
void print_POS_LST(POS_LST pos_lst);
void print_COUNT_LST(COUNT_LST cnt_lst);


Quantile quantile_ini(COUNT_LST cnt_lst, double p);

void quantile_update(LST lst,Quantile q,int pos,int read_cnt);
/* update the quantile given the new position <pos> and the number of reads <read_cnt>  at this position
 * <q> is the old quantile which will be updated.
 * this function will update the list <lst>.
 * */

void singularity_rm(int *read_pos, int num_pos, int win_size, double p, double multple); 
/* The major function of this file
 * read_pos: num_pos by 2 matrix, with its first column the positions of the read and its second column the number of reads at this position. In other words, each row of read_pos is a bin of 1 bp.
	     this should be ordered increasingly according to the genomic positions
 * num_pos: number of positions that has at least one read
 * win_size: the window size used to calculate the quantile
 * p: the pth quantile is used
 * multple: the genomic position whose read count will be set as <pth quantile>*multple if its read count is larger than that number.
 * */


#endif

