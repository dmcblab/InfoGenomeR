/*written by Ruibin Xi April 26 2010*/

#include "pos_cnt_lst.h"

static COUNT_LST_NODE cnt_node_create(int read_count);
/*creat a new node of type COUNT_LST_NODE given a read_count at certain position*/
static POS_LST_NODE pos_node_create(int pos, int read_cnt);

static void LST_add(LST lst, int pos, int read_cnt);
static void outlier_replace(int *bin, int i, Quantile q, double multpl);
/*Determine if the ith position is an outlier, if it is, set its read count as the quantile*multpl+1 */

static COUNT_LST_NODE cnt_node_create(int read_count){
		COUNT_LST_NODE node = malloc(sizeof(struct COUNT_LST_NODE_t));
		node->num_reads = read_count;
                node->num_pos = 1;
                node->pre = NULL;
                node->next = NULL;
 		return node;
	}


static POS_LST_NODE pos_node_create(int pos, int read_cnt){
		POS_LST_NODE node = malloc(sizeof(struct POS_LST_NODE_t));
		node->pos = pos;
		node->read_cnt = read_cnt;
		node->pre = NULL;
		node->next = NULL;
		node->cnt_node = NULL;
		return node;
}


COUNT_LST COUNT_LST_create(int win_size){ 
		COUNT_LST ll = malloc(sizeof(struct COUNT_LST_t));
		ll -> w = win_size;
		ll -> head = NULL;
		ll -> tail = NULL;
		return ll;
	}


COUNT_LST_NODE  insert_CNT(COUNT_LST ll, int read_count){
		COUNT_LST_NODE node, node_cur;
		
		if(ll==NULL) 
			{ printf("Error in \"insert_CNT\": the list is not initialized.\n"); 
			  exit(1);
			}
		
		if(ll->head==NULL)
			{ node = cnt_node_create(read_count); /*The first node in the list*/
			  ll->head = ll->tail = node;
			  return node;
			}

		node_cur = ll->head;
		while(node_cur!=NULL&&node_cur->num_reads > read_count)
			{ node_cur = node_cur->next;
			}/*find the right position to insert a node or to update the identified node*/

		if(node_cur==NULL) /*the new read_count is the smallest*/
			{ node = cnt_node_create(read_count);
			  node->pre = ll->tail;
			  ll->tail->next = node;
			  ll->tail = node;
			}
		else if(node_cur!=NULL&&node_cur->num_reads == read_count)
			{ node_cur->num_pos++;
			  node = node_cur;
			}
		else
			{  node = cnt_node_create(read_count); /*there is no such node whose num_reads is read_count, create a new node*/
			  if(node_cur==ll->head)
				{ ll->head->pre = node;
				  node->next = ll->head;
				  ll->head = node;
				}
			  else 
				{ node_cur->pre->next = node;
				  node->pre = node_cur -> pre;
				  node_cur->pre = node;
				  node->next = node_cur;
				}
			}

		 return node;
	}

void delete_CNT(COUNT_LST ll, COUNT_LST_NODE node)
	{ 
	  if(ll==NULL)
               { printf("Error in \"delete_CNT\": the list is not initialized.\n");
                 exit(1);
               }

	  node->num_pos--;

	  if(node->num_pos<=0)
		{ if(ll->head==node&&ll->tail==node) {ll->head=NULL; ll->tail=NULL;}
		  if(ll->head==node&&ll->tail!=node) 
			{ ll->head = node->next; node->next->pre = NULL;}
		  else if(ll->head!=node&&ll->tail==node) 
			{ ll->tail = node->pre;  node->pre->next = NULL;}
		  else
			{ node->pre->next = node->next;
			  node->next->pre = node->pre;
			}
		  free(node);
		}

	  return;
	}


POS_LST POS_LST_create()
	{ POS_LST ll = malloc(sizeof(struct POS_LST_t));
	  ll->head = NULL;
	  ll->tail = NULL;
	  ll->len = 0;
	  return ll;
	}

void push_POS(POS_LST ll, int pos, int read_cnt)
	{ POS_LST_NODE node;
	  
	  if(ll==NULL)
		{ printf("Error in \"push_POS\": list not initialized.\n");exit(1);}	
	  node = pos_node_create(pos,read_cnt);

	  if(ll->head==NULL)
		{ ll->head = node;
		  ll->tail = node;
		}
	  else {
	  	  node -> next = ll->head;
	  	  ll->head->pre = node; 
	  	  ll->head = node;
		}
	  ll->len++;
	  return;
	}

void del_tail_POS(POS_LST ll)
	{ POS_LST_NODE node;

          if(ll==NULL)
                { printf("Error in \"del_tail_POS\": list not initialized.\n");exit(1);}
	  
	  node = ll->tail;
	  if(node == NULL) return;

	  if(node == ll->head) {ll->head = NULL;ll->tail=NULL;}
	  else 
		{ ll->tail = node->pre;
		  ll->tail->next = NULL;
		}
	  ll->len--;
	  free(node);

	  return;
	}








LST LST_create(int w)
        { LST ll = malloc(sizeof(struct LST_t));
          ll->cnt_lst = COUNT_LST_create(w);
          ll->pos_lst = POS_LST_create();
          return ll;
        }


static void LST_add(LST lst, int pos, int read_cnt)
	{ COUNT_LST_NODE node;

	  push_POS(lst->pos_lst,pos,read_cnt);
	  node = insert_CNT(lst->cnt_lst,read_cnt);
	  lst->pos_lst->head->cnt_node = node;	  
	}


void lst_initialize(int *bin, int nbins, LST lst, int win_size)
        { int i;

          if(lst==NULL)
                { printf("Error in \"lst_initialize\": one of the linked list not initialized.\n");
                  exit(1);
                }
          lst->cnt_lst->w = win_size;

          i = 0;
          while(i<nbins&&i<win_size)
                { LST_add(lst,bin[2*i],bin[2*i+1]);
                  i++;
                }
          if(i<win_size) lst->cnt_lst->w = i;

          return;
        }


void print_POS_LST(POS_LST pos_lst)
        { POS_LST_NODE node;

          node = pos_lst->tail;


          while(node!=NULL)
                { printf("%d\t%d\t%d\n",node->pos,node->read_cnt,node->cnt_node->num_reads);
                  node = node->pre;
                }
          return;
        }


void print_COUNT_LST(COUNT_LST cnt_lst)
	{ COUNT_LST_NODE node;

	  node = cnt_lst -> head;

	  while(node!=NULL)
		{ printf("%d\t%d\n",node->num_reads,node->num_pos);
		  node = node->next;
		}
	  return;
	}


Quantile quantile_ini(COUNT_LST cnt_lst, double p)
	{ double q;/*q = w*(1-p), w is the window size*/
	  int k;
	  Quantile qntl = malloc(sizeof(struct Quantile_t));

	  if(cnt_lst==NULL)
		{ printf("Error in \"quantile_ini\": COUNT_LST not initialized.\n");
		  exit(1);
		}

	  qntl->p = p;
	  qntl->w = cnt_lst->w;
	  q = ((double) cnt_lst->w)*p;
	  k = 0;
	  qntl->cnt_node = cnt_lst->tail;  
	  while(((double)qntl->cnt_node->num_pos+k)-q<-1e-20 && qntl->cnt_node->pre!=NULL)
		{ k += qntl->cnt_node->num_pos;
		  qntl->cnt_node = qntl->cnt_node->pre; 
		}
	  qntl->k = k;
	  return qntl;
	}


void quantile_update(LST lst,Quantile q,int pos,int read_cnt)
	{ int diff; 
	  COUNT_LST_NODE head, tail;
	  head = lst->cnt_lst->head;
	  tail = lst->cnt_lst->tail;

	  //printf("\nold read count = %d\t new read count = %d\n",lst->pos_lst->tail->read_cnt,read_cnt);
	  //printf("old quantile: %d\n",q->cnt_node->num_reads);
	  //printf("before deletion: q->k = %d\n",q->k);
	  /*delete the old genomic position (slide the window by 1 to the right, first step)*/
	  if(lst->pos_lst->tail->read_cnt == q->cnt_node->num_reads)
		{ if(q->cnt_node->num_pos == 1)
			{ if(head == tail)
				{ q->cnt_node=NULL;q->k=0;}
			  else if(q->cnt_node==head)
				{ q->cnt_node = q->cnt_node->next; q->k -= q->cnt_node->num_pos;}
			  else {q->cnt_node = q->cnt_node->pre;}
			  //printf("yes, there is only one in this class\n");
			}
		  delete_CNT(lst->cnt_lst,lst->pos_lst->tail->cnt_node);
		  del_tail_POS(lst->pos_lst);
		  //printf("Yes, equal\n");
		}
	  else 	{ diff = (lst->pos_lst->tail->read_cnt >= q->cnt_node->num_reads)? 0:-1;
		  q->k += diff;

		  delete_CNT(lst->cnt_lst,lst->pos_lst->tail->cnt_node);
                  del_tail_POS(lst->pos_lst);
		  //printf("No, unequal, diff=%d\n",diff);
		}
	  //printf("updated quantile: %d\n",q->cnt_node->num_reads);
	  //printf("before insertion: q->k = %d\n",q->k);


	  /*insert the new genomic postion (slide the window by 1 to the right, second step)*/
	  LST_add(lst,pos,read_cnt);

	  //printf("q->w*q->p is %g\n",q->w*q->p);
	  if(q->cnt_node==NULL)
		{ q->cnt_node = lst->cnt_lst->head;}
	  else
		{ diff = (read_cnt >= q->cnt_node->num_reads) ? 0 : 1;
		  q->k += diff;
		  if(((double)q->k+q->cnt_node->num_pos) - q->w*q->p < -1e-20) /*In case of q->k+q->cnt_node->num_pos < q->w*q->p*/
			{ q->k += q->cnt_node->num_pos;
			  q->cnt_node = q->cnt_node->pre; 
			}
		  else if(((double)q->k) - q->w*q->p > -1e-20 && q->k!=0) /*In case of q->k >= q->w*q->p*/
			{ q->cnt_node = q->cnt_node->next;
			  q->k -= q->cnt_node->num_pos;
			}
		}
	  //printf("new quantile: %d\n",q->cnt_node->num_reads);
	  //printf("after insertion: q->k = %d\n\n",q->k);

	  return;
	}





static void outlier_replace(int *bin, int i, Quantile q, double multpl)
        { if(bin[2*i+1]>(q->cnt_node->num_reads)*multpl)
                { bin[2*i+1] = (int) (q->cnt_node->num_reads)*multpl;
                }
          return;
        }


void singularity_rm(int *bin, int nbins, int win_size, double p, double multple)
	{ int i,j,w_left,w_right;
	  LST lst;
          Quantile q;

	  if(p<0.0)
		{ printf("The probability specified for quantile calculation must be positive. (prob=%g)\n",p);
		  exit(1);
		}
	  if(multple<1.0)
		{ printf("the mutiplicity of the quantile MUST be larger than 1( m = %g)\n",multple);
		  exit(1);
		}
	  if(win_size<1)
		{ printf("The window size (=%d) must be larger than or equal to 1.\n",win_size);
		  exit(1);
		}
	  if(win_size > nbins)
		{ printf("Warning: the window size (=%d) is larger than the number of genomic positions (=%d) that have at least one reads.\n",win_size,nbins);
		}	

	  lst = LST_create(win_size);
	  lst_initialize(bin,nbins,lst,win_size);
	  q = quantile_ini(lst->cnt_lst,p);
	  w_left = (q->w-1)/2; /*the left window size*/
	  w_right = q->w - w_left -1; /*the right window size*/

	  i = 0;/*the central position in the window (except on the boundary of the chromosome)*/
	  outlier_replace(bin,i,q,multple);
	  for(i=1;i<nbins;i++)
		{ if(i+1>w_left&&i+1+w_right<nbins)
			{ j = i + (q->w-1)/2+1;
			  quantile_update(lst,q,bin[2*j],bin[2*j+1]);
			}
		  if(q->w!=lst->pos_lst->len) {printf("Error!\n");exit(1);}
		  outlier_replace(bin,i,q,multple);
		}
	  
	  return;
	}
