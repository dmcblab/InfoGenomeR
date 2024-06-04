#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

void
dump_line(FILE *in, FILE *out)
{
	int n;
	char c;

	n = 0;
	while (!feof(in)) {
		c = fgetc(in);
		if (c == '\n' || c == '\r' || c == EOF)
			break;
		++n;
		fprintf(out, "%c", c);
	}

	if (n > 0)
		fprintf(out, "\n");
}


void read_line(FILE *in)
{
        int n;
        char c;

        n = 0;
        while (!feof(in)) {
                c = fgetc(in);
                if (c == '\n' || c == '\r' || c == EOF)
                        break;
                ++n;
        }

}

int
main(int argc, char **argv)
{
	double p;
	FILE *f;
	FILE *f1;
	struct timeval tv;
	long seed;

	if (argc != 2){
		fprintf(stderr, "usage: %s <probability that line goes to output file>\n",argv[0]);
		fprintf(stderr, "     the input should be <stdin> and the output is <stdout>\n");
		exit(1);
		}

	p = strtod(argv[1], NULL);

	f = stdin;

	//f1 = fopen(argv[2], "w");
	//if (!f1) err(1, "main: fopen");
	f1 = stdout;

	gettimeofday(&tv, NULL);
	seed = tv.tv_sec * 1000000 + tv.tv_usec;
	srand48(seed);

	fprintf(stderr,"p=%lf,seed=%ld\n",p, seed);

	while (!feof(f)) {
		if (drand48() < p)
			dump_line(f, f1);
		else read_line(f);
	}

	//fclose(f1);

	return 0;
}
