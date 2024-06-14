#ifndef GetUniqueBamRead
#define GetUniqueBamRead

#define _GNU_SOURCE     /* Expose declaration of tdestroy() */
#include <search.h>

#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "../bam.h"
#include "../sam.h"



int get_Unique_args(char *args);

void report_unique_read(bam_header_t *in,bam1_t *b);
void summarizeReadLength(bam1_t *b);
void reportReadLength(FILE *out);

#endif

