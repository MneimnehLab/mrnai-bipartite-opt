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
// FILE * fOutput;
int wType = 1;
int gap_size;
extern double ***** rnaupCollections;

extern int * blahblahptr();
extern CmdLineArgs * parseArgs(int argc, char *argv[]);
extern int readAndTokenizeRNAs(char **strArray1);
extern void readEvenOdd();
extern void tokenizeNames(int num);
extern void doCachingStart();
extern void computeMatchingMatrix(int numOfRNA);
extern void computeMatchingMatrixFromRNAup(int numOfRNA, char * rnaupOut, int GU);
extern void computeMatchingMatrixFromFile(int numOfRNA, char * fileName, char * rnaupOut);
extern Pointer * makeSampleBinStructJustTwo();
extern Pointer * makeSomeBinStruct(OrigRNASeq ** origRNASequences, int numOfRNA);
extern void printChain(Pointer * head);
extern void determineStruct(Pointer * head);
extern int getBBHeight(Pointer * head);
extern void printFlatStruct2File(Pointer * head, FILE * file);
extern void printFlatInterStructure2FileGaps(Ensemble * ensemble, FILE * fp);
extern void printFlatInterStructure2File(Ensemble * ensemble, FILE * fp);
extern Ensemble * doSplitting(int levels, int k, double * energy, Pointer * head);
extern void drawThing(Ensemble * ensemble, OrigRNASeq ** origRNASequences);
extern void destroyStruct(Pointer * head);


void printMilliSecs()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);

	long t2 = ((tv.tv_usec + 1000000 * tv.tv_sec) - (lastTV.tv_usec + 1000000 * lastTV.tv_sec) ) / 1000;
	gettimeofday(&lastTV, NULL);	
	printf("\033[1;31mElapsed Time \033[0m = %ld \n", t2);
}

void printTime()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);

	long t2 = ((tv.tv_usec + 1000000 * tv.tv_sec) ) / 1000;
	printf("\033[1;31m !! Time !! \033[0m = %ld \n", t2);
}




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
	// printMilliSecs();

	// printf("hello\n");
	
	// fOutput = fopen(PROG_TYPE "results.txt", "w");
	// fprintf (fOutput, "Process started at: %s \n\n", asctime (timeinfo) );
	// fflush(fOutput);
	
	int i;
	for(i=0;i<10;i++)
		origRNASequences[i] = (OrigRNASeq *) malloc(sizeof(OrigRNASeq));

	//
	// /** 1. Read RNAs **/
	// numOfRNA = readAndTokenizeRNAs(strArray1);
	
	// /** 2. Get the (even|odd) type of each RNA **/
	// tokenizeNames(numOfRNA);
	
	
	// readEvenOdd();

	numOfRNA = getInput(strArray1);

	printf("\n");
	for(i=0;i<numOfRNA;i++)
	{
		printf("%d name:     %s\n", i, origRNASequences[i]->name);
		printf("%d sequence: %s\n", i, origRNASequences[i]->string);
		printf("%d length:   %d\n", i, origRNASequences[i]->originalLength);
	}


	
	// doCachingStart();

	//computeMatchingMatrix(numOfRNA, args->rnaupOut);
	
	// if(strcmp(args->fileName, "default") == 0)
	computeMatchingMatrixFromRNAup(numOfRNA, args->rnaupOut, args->GU);
	// else
	// computeMatchingMatrixFromFile(numOfRNA, args->fileName, args->rnaupOut);
	
	printf("Have computed matching matrix\n");


	return 0;


	
	printTime();
	
	/* Print given regions energy 
	printf("\nRegion values::\n");
	printf("%f\n", rnaupCollections[0][58][38][8][8]);
	printf("%f\n", rnaupCollections[0][14][15][3][3]);
	printf("%f\n\n", rnaupCollections[0][4][4][3][3]);
	
	printf("Region values::\n");
	printf("%f\n", rnaupCollections[0][58][38][8][8]);
	printf("%f\n", rnaupCollections[0][17][20][2][2]);
	printf("%f\n", rnaupCollections[0][14][15][3][3]);
	printf("%f\n\n", rnaupCollections[0][4][4][3][3]);
	printf("\n%f\n\n", rnaupCollections[0][17][20][6][8]);
	*/
	
	
	FILE * randOut = fopen(PROG_TYPE "randout.txt", "w");
	double fe = 0;
	int zx;
		
	
	FILE * file = fopen(PROG_TYPE "all_perms_flat.txt", "r");
	FILE * fw = fopen(PROG_TYPE "all_perms_flat_output.dat", "w");
	int lineC = 0;
	
	//Pointer * head = (Pointer *) makeSampleBinStruct4(origRNASequences, numOfRNA);
	//Pointer * head = (Pointer *) makeSampleBinStruct5(origRNASequences, numOfRNA);
	
	
	
	//Pointer * head = (Pointer *) makeSampleBinStruct15(origRNASequences, numOfRNA);
	//Pointer * head = (Pointer *) testConvgOutStruct1(origRNASequences, numOfRNA);
	
//Pointer * head = (Pointer *) makeRandomBinStruct(origRNASequences, numOfRNA);
	
	//Pointer * head = (Pointer *) makeSampleBinStruct5(origRNASequences, numOfRNA);
	
	//Pointer * head = (Pointer *) makeSampleBinStructJustTwo();
	//Pointer * head = (Pointer *) makeSampleBinStructJustTwoStartOdd();
	printf("num of rna = %d\n", numOfRNA);
	Pointer * head;
	if(numOfRNA == 2)
		head = (Pointer *) makeSampleBinStructJustTwo();
	//else if(numOfRNA == 3)
	//	head = (Pointer *) makeU6_U4U2();
	else
		//head = (Pointer *) justFourOddToEven();
		head = (Pointer *) makeSomeBinStruct(origRNASequences, numOfRNA);
		
	//Pointer * head = (Pointer *) justFourOddToEven2_31_0();
	
//Pointer * head = (Pointer *) somethingWrongInMixedCnvg_1();
//Pointer * head = (Pointer *) testMixedConvg1();
//Pointer * head = (Pointer *) somethingWrongInCnvg2();
	//Pointer * head = (Pointer *) makeBinStruct_Highest();
	//Pointer * head = (Pointer *) somethingWrongInCnvg_12();	// = (Pointer *) createFromString("72145630");
	//Pointer * head = (Pointer *) somethingWrongInMixedCnvg_2();	
	
	//Pointer * head = (Pointer *) makeBinStruct_Highest();	
	//Pointer * head = (Pointer *) somethingWrongInCnvg_12();	
	
	//Pointer * head = (Pointer *) makeSampleBinStructCopACopT();
	
	
	//Pointer * head = (Pointer *) nice6_2();
	
	printChain(head);
	determineStruct(head);
	
	int height = getBBHeight(head);
	printf("Height = %d \n", height);


	fprintf(randOut, "Starter: ");
	printFlatStruct2File(head, randOut);
	printFlatStruct2File(head, fw);
	
	fprintf(randOut, "\n");
	fflush(randOut);
	
	/* 3. Somewhere here, we will put the tree shuffling stuff */
	
	totalLevels = height;
	
	/*** This part is specific to a confguration of the tree structure ***/
	/*** In the future, it will be done automatically ***/
	
	
	/*** End - This part is specific to a confguration of the tree structure ***/
	
	
	/** At this point we know the number of levels and the layout of the RNAs. 
		We should now run the following code for subsets of levels (of size k) **/

	
	

	int k = args->k;

	if(height < k) k = height;

	// fprintf(fOutput, "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n");
	// fprintf(fOutput, "Working with struct 1: \n");
	// printFlatStruct2File(head, fOutput);
	// fprintf(fOutput, "\n\n");
	// fprintf(fOutput, "Starting split with k = %d\n\n", k);
	// fflush(fOutput);
	
	
	double energy;
	Ensemble * ensemble = (Ensemble * ) doSplitting(totalLevels, k, &energy, head);

	fprintf(randOut, "\n");
	fprintf(randOut, "Energy: %f \t| First = %d\n\n", energy, ensemble->first);
	
	//fprintf(fw, "\tEnergy: %f\tFirst = %d\t", energy, ensemble->first);
	fprintf(fw, "\t%f\t%d\t", energy, ensemble->first);
	printFlatInterStructure2File(ensemble, fw);
	fprintf(fw, "\n");

	drawThing(ensemble, origRNASequences);
	
	fflush(fw);
	
	destroyStruct(head);
	
	fprintf(randOut, "\n");
	fprintf(randOut, "\n");

	fflush(randOut);	


	fclose ( file );
	fclose ( fw );
	
	// fclose(fOutput);
	fclose(randOut);
	free(args);
	
	
	printTime();
	
	
	return 0;
}



