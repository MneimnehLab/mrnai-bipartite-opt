#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <float.h> 
#include <stdarg.h>
#include <fcntl.h> 
#include <time.h>
#include <sys/time.h>
#include "extractor/rnaseq.h"
#include "extractor/bins.h"
#include "extractor/datastruct.h"
#include "extractor/bins_flat.h"

#define PROG_TYPE "findone_"

long startT;
long lastT = 0;
struct timeval lastTV;

double INF1 = 9998;

int numOfRNA = 0;
char ** strArray1;
int setup[10][10];
int parentOf[10][10];
int totalLevels = 0;
int parallel = 0;
int fillGaps = 1;

OrigRNASeq * origRNASequences[10];
int wType = 1;
int gap_size;

extern CmdLineArgs * parseArgs(int argc, char *argv[]);
extern void computeMatchingMatrixFromRNAup(int numOfRNA, char * rnaupOut, int GU);

int main(int argc, char *argv[])
{
	
	CmdLineArgs * args;
	args = (CmdLineArgs * ) parseArgs(argc, argv);
	fillGaps = args->fillGaps;
	wType = args->winSize;
	gap_size = args->gapSize;	

	/* Initialize */
	time_t start,end, curr;
	double timeDiff;
	time (&start);
	struct tm * timeinfo;
	timeinfo = localtime ( &start );

	gettimeofday(&lastTV, NULL);

	int i;
	for(i=0;i<10;i++)
		origRNASequences[i] = (OrigRNASeq *) malloc(sizeof(OrigRNASeq));

	numOfRNA = getInput(strArray1);

	printf("\n");
	for(i=0;i<numOfRNA;i++)
	{
		printf("%d name:     %s\n", i, origRNASequences[i]->name);
		printf("%d sequence: %s\n", i, origRNASequences[i]->string);
		printf("%d length:   %d\n", i, origRNASequences[i]->originalLength);
	}

	computeMatchingMatrixFromRNAup(numOfRNA, args->rnaupOut, args->GU);
	printf("Have computed matching matrix\n");

	return 0;
}



