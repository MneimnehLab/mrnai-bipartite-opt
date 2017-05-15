#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <float.h> 
#include <stdarg.h>
#include <sys/time.h>
#include "rnaseq.h"

extern OrigRNASeq * origRNASequences[10];
extern int wType;
int matchingMatrix[10][10];
double ***** rnaupCollections;

// extern void doRNAupAndReverse(char * evenString, char * oddString, int location, double ***** rnaupCollections, char * rnaupOut1, int GU);

void doSubAddFor2RNAs_twoW(int loc, int n1, int n2, int rna1, int rna2);
void doSubAddFor2RNAs_sameW(int loc, int n1, int n2, int rna1, int rna2);
void readPirna(int evenId, int oddId, int n1, int n2, int location, double ***** rnaupCollections, char * fileName2);
void doSubAddFor2RNAs_sameW(int loc, int n1, int n2, int rna1, int rna2);


//For legacy reasons
void computeMatchingMatrix(int numOfRNA)
{
}
void computeMatchingMatrixFromRNAup(int numOfRNA, char * rnaupOut, int GU)
{
	printf("Computing matrix from RNAup.... \n\n");
	int xi;
	// for(xi=0;xi<numOfRNA;xi++)
	// 	printf("Z1 %d: %s \n", xi, origRNASequences[xi]->name);


	struct timeval startTV;
	gettimeofday(&startTV, NULL);

	int i;
	short evenSeq[10], e = 0;
	short oddSeq[10], o = 0;

	// for(xi=0;xi<numOfRNA;xi++)
	// 	printf("Z2 %d: %s \n", xi, origRNASequences[xi]->name);



	//initialize matrix
	for(i=0;i<10*10;i++)
		// matchingMatrix[i/10][i%10] = -9;	
		matchingMatrix[0][i] = -9;	

	//find even or odd
	for(i=0;i<numOfRNA;i++)
	{
		if(origRNASequences[i]->type == 0)
		{
			// printf("Even %d = %d \n", e, i);
			evenSeq[e++] = i;
		}
		else if(origRNASequences[i]->type == 1)
		{
			// printf("Odd %d = %d \n", o, i);
			oddSeq[o++] = i;	
		}
	}
	

	//do RNAup for each even/odd pair
	int oi,ei, pairs = o * e;
	int maxEO = e > o ? e : o;

	double ***** rnaupCollections2;

	rnaupCollections = (double*****)malloc(sizeof(double ****) * pairs);
	rnaupCollections2 = rnaupCollections;

	i = 0;
	for(oi=0;oi<o;oi++)
	{
		for(ei=0;ei<e;ei++)
		{

			int even = evenSeq[ei];
			int odd  = oddSeq[oi];

			printf("Pair even: %d and odd: %d, num = %d \n", even, odd, i);

			matchingMatrix[even][odd] = i;
			// matchingMatrix[odd][even] = i;


			printf("Even name = %s \n", origRNASequences[even]->name);
			printf("Odd name = %s \n", origRNASequences[odd]->name);
			// printf("Even str = %s \n", origRNASequences[even]->string);
			// printf("Odd str = %s \n", origRNASequences[odd]->string);


			char rnaupOutExtend[40];

			sprintf (rnaupOutExtend, "%s-%d_%d", rnaupOut, even, odd);
			doRNAupAndReverse(origRNASequences[even]->string, origRNASequences[odd]->string, i, rnaupCollections, rnaupOutExtend, GU);
			
			// printf("printing... \n");
			
			int x,y;
			char fileName[40];
			// sprintf (fileName, "%s-file-%d_%d.dat", rnaupOut, even, odd);
			// FILE * fp = fopen(fileName, "w");
			
			int evenLen = strlen(origRNASequences[even]->string);
			int oddLen = strlen(origRNASequences[odd]->string);
			
			if(wType == 1)
				doSubAddFor2RNAs_sameW(i, evenLen, oddLen, even, odd);
			else
				doSubAddFor2RNAs_twoW(i, evenLen, oddLen, even, odd);
			
			i++;


		}
	}

	int j;
	printf("Matching matrix: \n");
	for(i=0;i<e+o;i++)
	{
		for(j=0;j<e+o;j++)
		{
			printf("%2d   ", matchingMatrix[i][j]);
			//rnaupCollections
			int rel = matchingMatrix[i][j];
			if(rel != -9)
			{	


			}
		}
		printf("\n");
	}


	printf("\nFileID -> Even Odd:\n");
	for(i=0;i<e+o;i++)
	{
		for(j=0;j<e+o;j++)
		{
			if(matchingMatrix[i][j] != -9)
				printf("%d -> %d %d\n", matchingMatrix[i][j], i, j);
		}
	}

	struct timeval endTV;	
	gettimeofday(&endTV, NULL);		
	long t2 = ((endTV.tv_usec + 1000000 * endTV.tv_sec) - (startTV.tv_usec + 1000000 * startTV.tv_sec) ) / 1000;

	printf("\033[0;35mTotal Time for %d RNAups \033[0m = %ld \n", pairs, t2);

}



void computeMatchingMatrixFromFile(int numOfRNA, char * fileName, char * rnaupOut)
{
	printf("Read weights from file.... \n\n");
	int xi;
	for(xi=0;xi<numOfRNA;xi++)
		printf("Z1 %d: %s \n", xi, origRNASequences[xi]->name);


	struct timeval startTV;
	gettimeofday(&startTV, NULL);

	int i;
	short evenSeq[10], e = 0;
	short oddSeq[10], o = 0;

	for(xi=0;xi<numOfRNA;xi++)
		printf("Z2 %d: %s \n", xi, origRNASequences[xi]->name);


	//initialize matrix
	for(i=0;i<100;i++)
		matchingMatrix[i/10][i%10] = -9;	

	//find even or odd
	for(i=0;i<numOfRNA;i++)
	{
		if(origRNASequences[i]->type == 0)
		{
			printf("Even %d = %d \n", e, i);
			evenSeq[e++] = i;
		}
		else if(origRNASequences[i]->type == 1)
		{
			printf("Odd %d = %d \n", o, i);
			oddSeq[o++] = i;	
		}
	}
	for(xi=0;xi<numOfRNA;xi++)
		printf("Z3 %d: %s \n", xi, origRNASequences[xi]->name);

	//do RNAup for each even/odd pair
	int oi,ei, pairs = o * e;

	double ***** rnaupCollections2;

	rnaupCollections = (double*****)malloc(sizeof(double ****) * pairs);
	rnaupCollections2 = rnaupCollections;

	i = 0;
	for(oi=0;oi<o;oi++)
	{
		for(ei=0;ei<e;ei++)
		{

			int even = evenSeq[ei];
			int odd  = oddSeq[oi];

			printf("Pair even: %d and odd: %d, num = %d \n", even, odd, i);

			matchingMatrix[even][odd] = i;
			matchingMatrix[odd][even] = i;


			printf("Even name = %s \n", origRNASequences[even]->name);
			printf("Odd name = %s \n", origRNASequences[odd]->name);
			int evenLen = strlen(origRNASequences[even]->string);
			int oddLen = strlen(origRNASequences[odd]->string);
			
			
			readPirna(origRNASequences[even]->origId, origRNASequences[odd]->origId, evenLen, oddLen, i, rnaupCollections, fileName);
			
			if(wType == 1)
				doSubAddFor2RNAs_sameW(i, evenLen, oddLen, even, odd);
			else
				doSubAddFor2RNAs_twoW(i, evenLen, oddLen, even, odd);
			
			
			
			/*
			printf("printing... \n");
			
			int x,y;
			char fileName[20];
			sprintf (fileName, "file-%d_%d.dat", even, odd);
			
			FILE * fp = fopen(fileName, "w");
			
			
			
			
			
			for(x=0;x<=evenLen;x++)
			{
				for(y=0;y<=oddLen;y++)
				{
					int w, w2;
					for(w=1;w<=25;w++)
					{
						for(w2=1;w2<=25;w2++)
						{
							fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%f \n", i, x, y, w, w2, rnaupCollections[i][x][y][w][w2]);
						}
					}
				}
			}
			
			printf("printed \n");
			fclose(fp);
			*/
			
			
/*
			int x,y,w;
			for(x=1;x<50;x++)
				for(y=1;y<50;y++)
					for(w=1;w<=25;w++)
						printf("rnaupCollections[%d][%d][%d][%d] = %f \n", i, x,y,w, rnaupCollections[i][x][y][w]);
*/
			i++;


		}
	}

	int j;
	printf("Matching matrix: \n");
	for(i=0;i<10;i++)
	{
		for(j=0;j<10;j++)
		{
			printf("%d   ", matchingMatrix[i][j]);
			//rnaupCollections
			int rel = matchingMatrix[i][j];
			if(rel != -9)
			{	/*
				int x,y;
				char fileName[20];
				sprintf (fileName, "d:\\file-%d_%d.dat", x, y);
				FILE * fp = fopen(fileName, "w");
				
				for(x=0;x<MAX_RNA_SIZE;x++)
				{
					for(y=0;y<MAX_RNA_SIZE;y++)
					{
						int w;
						for(w=1;w<=25;w++)
						{
							fprintf(fp, "%f \n", rnaupCollections[rel][x][y][w]);
						}
					}
				}*/
			}
		}
		printf("\n");
	}

	struct timeval endTV;	
	gettimeofday(&endTV, NULL);		
	long t2 = ((endTV.tv_usec + 1000000 * endTV.tv_sec) - (startTV.tv_usec + 1000000 * startTV.tv_sec) ) / 1000;

	printf("\033[0;35mTotal Time for %d RNAups \033[0m = %ld \n", pairs, t2);

}



void doSubAddFor2RNAs_twoW(int loc, int n1, int n2, int rna1, int rna2)
{
	// printf("In subadditive for two W... \n");

	int WMAX = 25;
	int w, i, j, w2;
	
	for(w=1;w<=WMAX;w++)
	{
		for(w2=1;w2<=WMAX;w2++)
		{
		
			//FOR a specific window size, do:
			
			
			for(i=0;i<n1;i++)
			{
				for(j=0;j<n2;j++)
				{
					//printf("Here 1\n");
					double windowEnergy = rnaupCollections[loc][i][j][w][w2];
					//printf("Here 2\n");
					if(windowEnergy >= 0) 
						continue;
					//printf("Here 3\n");
					//define end points of the RNAs
					int rna1_l = i-w;
					int rna2_l = j-w2;
					int rna1_r = i;
					int rna2_r = j;
					//printf("Here 4\n");
					if(rna1_l < 0 || rna2_l < 0) 
						continue;
					
					//printf("Here 5\n");
					//find a split point
					int p = 1, q;
					for(p=1;p<w;p++)
					{

						for(q=1;q<w2;q++)
						{
							//printf("Here 6\n");
							//printf("w:%d, w2:%d, i:%d, j:%d, p:%d, q:%d \n", w,w2,i,j,p,q);
							//printf("rna1_l:%d, rna2_l:%d, rna1_r:%d, rna2_r:%d\n",  rna1_l, rna2_l, rna1_r, rna2_r);
							//printf("rna1_l+p:%d, rna2_l:%d, rna1_r:%d, rna2_r:%d\n",  rna1_l+p, rna2_l, rna1_r, rna2_r);
						
							//left half:    [rna1_l, rna1_l+p] : [rna2_l, rna2_l+q];;; win size = p
							//right half: [rna1_l+p+1, rna1_r] : [rna2_l+q+1, rna2_r];;; win size = (rna1_r) - (rna1_l+p+1) = rna1_r-rna1_l-p-1
							
							double leftChild  = rnaupCollections[loc][rna1_l+p][rna2_l+q][p][q];
							//printf("A \n");
							double rightChild = rnaupCollections[loc][rna1_r][rna2_r][rna1_r-rna1_l-p-1][rna2_r-rna2_l-q-1]; //win = w - p - 1;
							//printf("B \n");
							
							if(leftChild > 0 || rightChild > 0)
								continue;
							
							//printf("C \n");
							double sumOfChildren = leftChild + rightChild;
							//printf("D \n");
							
							if(windowEnergy >= sumOfChildren)	//Error!
							{
								//printf("E \n");
								
								// printf("Not sub-additive!\n");
								// printf("RNA1 : %d, [%d,%d] \n", rna1, rna1_l, rna1_r);
								// printf("RNA2 : %d, [%d,%d] \n", rna2, rna2_l, rna2_r);
								// printf("Energy: %f \n", windowEnergy);
								
								// printf("Left  Child [%d,%d] : [%d,%d], Energy = %f \n", rna1_l, rna1_l+p, rna2_l, rna2_l+p, leftChild);
								// printf("Right Child [%d,%d] : [%d,%d], Energy = %f \n", rna1_l+p+1, rna1_r, rna2_l+p+1, rna2_r, rightChild);
								// printf("\n");
								
								rnaupCollections[loc][rna1_l+p][rna2_l+p][p][q] = 999.0;
								rnaupCollections[loc][rna1_r][rna2_r][rna1_r-rna1_l-p-1][rna2_r-rna2_l-q-1] = 999.0; 
								//printf("F \n");
							}
							
						}
						
					}
					
				}
			}
			
		}
		
	}
	
	// printf("Finished subadditive...\n");
	
}



void doSubAddFor2RNAs_sameW(int loc, int n1, int n2, int rna1, int rna2)
{
	// printf("In subadditive for same W... \n");

	int WMAX = 25;
	int w, i, j, w2;
	
	for(w=1;w<=WMAX;w++)
	{
			//FOR a specific window size, do:
			
			
			for(i=0;i<n1;i++)
			{
				for(j=0;j<n2;j++)
				{
					double windowEnergy = rnaupCollections[loc][i][j][w][w];
					if(windowEnergy >= 0) 
						continue;
						
					//define end points of the RNAs
					int rna1_l = i-w;
					int rna2_l = j-w;
					int rna1_r = i;
					int rna2_r = j;
					
					if(rna1_l < 0 || rna2_l < 0) 
						continue;
					
					
					//find a split point
					int p = 1, q;
					for(p=1;p<w;p++)
					{
							q = p;
							//printf("w:%d, w2:%d, i:%d, j:%d, p:%d, q:%d \n", w,w2,i,j,p,q);
							//printf("rna1_l:%d, rna2_l:%d, rna1_r:%d, rna2_r:%d\n",  rna1_l, rna2_l, rna1_r, rna2_r);
							//printf("rna1_l+p:%d, rna2_l:%d, rna1_r:%d, rna2_r:%d\n",  rna1_l+p, rna2_l, rna1_r, rna2_r);
						
							//left half:    [rna1_l, rna1_l+p] : [rna2_l, rna2_l+q];;; win size = p
							//right half: [rna1_l+p+1, rna1_r] : [rna2_l+q+1, rna2_r];;; win size = (rna1_r) - (rna1_l+p+1) = rna1_r-rna1_l-p-1
							
							double leftChild  = rnaupCollections[loc][rna1_l+p][rna2_l+q][p][q];
							//printf("A \n");
							double rightChild = rnaupCollections[loc][rna1_r][rna2_r][rna1_r-rna1_l-p-1][rna2_r-rna2_l-q-1]; //win = w - p - 1;
							//printf("B \n");
							
							if(leftChild > 0 || rightChild > 0)
								continue;
							
							//printf("C \n");
							double sumOfChildren = leftChild + rightChild;
							//printf("D \n");
							
							if(windowEnergy >= sumOfChildren)	//Error!
							{
								//printf("E \n");
								
								// printf("Not sub-additive!\n");
								// printf("RNA1 : %d, [%d,%d] \n", rna1, rna1_l, rna1_r);
								// printf("RNA2 : %d, [%d,%d] \n", rna2, rna2_l, rna2_r);
								// printf("Energy: %f \n", windowEnergy);
								
								// printf("Left  Child [%d,%d] : [%d,%d], Energy = %f \n", rna1_l, rna1_l+p, rna2_l, rna2_l+p, leftChild);
								// printf("Right Child [%d,%d] : [%d,%d], Energy = %f \n", rna1_l+p+1, rna1_r, rna2_l+p+1, rna2_r, rightChild);
								// printf("\n");
								
								rnaupCollections[loc][rna1_l+p][rna2_l+p][p][q] = 999.0;
								rnaupCollections[loc][rna1_r][rna2_r][rna1_r-rna1_l-p-1][rna2_r-rna2_l-q-1] = 999.0; 
								//printf("F \n");
							}
							
					
						
					}
					
				
			}
			
		}
		
	}
	
	// printf("Finished subadditive...\n");
	
}

void readPirna(int evenId, int oddId, int n1, int n2, int location, double ***** rnaupCollections, char * fileName2)
{
	//allocate memory for rnaWins
	double **** rnaWins;	//rnaWins[1st rna][2nd rna][w1][w2];
	double INF_1 = 9999;
	rnaWins = (double ****) malloc(sizeof (double ***) * (n1 + 1));
	int i,j,k;
	for (i = 0; i <= n1; i++)
	{
		rnaWins[i] = (double ***) malloc(sizeof (double **) * (n2 + 1));
		for (j = 0; j <= n2; j++)
		{
			rnaWins[i][j] = (double **) malloc(sizeof (double *) * (26));
			for (k = 0; k <= 25; k++)
			{
				rnaWins[i][j][k] = (double *) malloc(sizeof (double) * (26));
				int k2 = 0;
				for (k2 = 0; k2 <= 25; k2++)
				{
					rnaWins[i][j][k][k2] = INF_1;
				}
			}
		}
	}
	
	char fileName[50];

	printf("filename2 = *%s* \n", fileName2);
	if(!strcmp(fileName2, "default"))
	{
		printf("here 1\n");
		printf("Reading default file... \n");
		sprintf(fileName, "pirna_out_%d_%d.dat", evenId, oddId);
		printf("File name = %s \n\n", fileName);
	}
	else
	{	printf("here 2 \n");
		strcpy(fileName, fileName2);
		printf("Reading file %s... \n", fileName);
	}
	
	
	int count = 0;
	FILE * file = fopen(fileName, "rt");
	if ( file != NULL )
	{
		char line [ 500 ]; 
		while ( fgets ( line, sizeof line, file ) != NULL ) 
		{		
			
			int i1,j1,i2,j2;
			double energy;
			//printf("%s\n", line);
			sscanf (line, "%d\t%d\t%d\t%d\t%lf\n", &i1,&j1,&i2,&j2,&energy);
			
			int w1 = j1-i1;
			int w2 = j2-i2;
			rnaWins[j1+1][j2+1][w1][w2] = energy;
			
			//printf("%d, %d, %d, %d = %f \n", j1+1,j2+1,w1,w2, energy);
			//printf("%d, %d, %d, %d = %f \n", i1,j1,i2,j2, energy);
			
		}
		fclose ( file );
		printf("Read file \n\n");
		rnaupCollections[location] = rnaWins;
	}
	else
	{
		printf("File read error! \n");
	}
	
}


