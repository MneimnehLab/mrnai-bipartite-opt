#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <float.h> 
#include <stdarg.h>
#include "rnaseq.h"

extern OrigRNASeq * origRNASequences[10];
extern int parallel;

int myTokenize(char *line, char **strArray1);

int readAndTokenizeRNAs(char **strArray1)
{
	int i;
	printf("Enter sequences: ");
	
	/// Use this to input ///
	char *line___ = NULL;
	size_t size;
	int lres = getline(&line___, &size, stdin);
	char line[500];
	strcpy(line, line___);
	line[strlen(line___)-1] = '\0';
	//////////////////////////
	
	strArray1 = (char **) malloc(sizeof(char *) * 10);
	int numOfRNA = myTokenize(line, strArray1);
	// printf("num = %d\n", numOfRNA);
	for(i=0;i<numOfRNA;i++)
	{
		origRNASequences[i]->originalLength = strlen(origRNASequences[i]->string);
		// printf("orig len = %d \n", origRNASequences[i]->originalLength);
	}	
	/*
	for(i=0;i<numOfRNA;i++)
		printf("seq %d: %s\n", i, origRNASequences[i]->string);
	*/
	// printf("Number of RNAs = %d \n", numOfRNA);
	return numOfRNA;
}


int myTokenize(char *line, char **strArray1)
{
	//assume atleast 1 & exists
	int numOfRNA;
	char * str = line;
	char * pch;
	
	pch = strtok (str,"&");
	int i = 0;
	
	for (; pch != NULL; i++, pch = strtok (NULL, "&"))
	{
		origRNASequences[i]->string = (char*) malloc(sizeof(char) *
											(strlen(pch)+1));
		strcpy(origRNASequences[i]->string, pch);
		origRNASequences[i]->origId = i;
	}
	
	return i;
}



void readEvenOdd()
{
	printf("Enter even|odd:");
	//char * line = get_line(stdin);
	
	/// Use this to input ///
	char *line___ = NULL;
	size_t size;
	int lres = getline(&line___, &size, stdin);
	char line[200];
	strcpy(line, line___);
	line[strlen(line___)-1] = '\0';
	//////////////////////////

	
	int len = strlen(line);
	
	int i;
	int readingType = 0;
	
	for(i=0;i<len;i++)
	{
		// printf("chr: %c\n", line[i]);
		switch(line[i])
		{
			case '|':
			{
				readingType = 1;
				break;
			}
			case ',':
			{
				break;
			}
			default:
			{
				origRNASequences[line[i]-48]->type = readingType;
				// printf("Read %d, gave type %d \n", line[i]-48, readingType);
				// printf("Id was %d \n", origRNASequences[line[i]-48]->origId);
				break;
			}
		}
		
	}
	// printf("done\n");
}















void tokenizeNames(int num)
{
	printf("Enter names seperated by \'&\' \n" );
	//char * line = get_line(stdin);
	
	/// Use this to input ///
	char *line___ = NULL;
	size_t size;
	int lres = getline(&line___, &size, stdin);
	char line[200];
	strcpy(line, line___);
	line[strlen(line___)-1] = '\0';
	//////////////////////////
	

	// assume atleast 1 & exists
	char * str = line;
	char * pch = strtok (str,"&");
	int i = 0;
	strcpy(origRNASequences[i]->name, pch);
	i++;
	
	while (i < num)
	{
		//printf ("%s\n",pch);
		pch = strtok (NULL, "&");
		strcpy(origRNASequences[i]->name, pch);
		i++;
	}
}






int blahblah()
{
	int x = 53;
	return x;
}

int * blahblahptr()
{
	int * p = (int*)malloc(sizeof(int));
	printf("p ptr = %p\n", p);
	*p = 65;
	return p;
}


CmdLineArgs * parseArgs(int argc, char ** argv)
{
	CmdLineArgs * args = (CmdLineArgs *)malloc(sizeof(CmdLineArgs));
	// printf("args ptr = %p\n", args);
	//default args values:
	args->parallel = 0;
	args->k = 2;
	args->trials = 1;
	args->fillGaps = 1;
	args->gapSize = 4;
	args->winSize = 1;
	args->GU = 1;
	strcpy(args->fileName, "default");
	strcpy(args->rnaupOut, "default");

	int i;
	int k, t, gbool, gsize, w, gu;
	for (i=1; i<argc; i++) 
	{
		if (argv[i][0]=='-') 
  			switch ( argv[i][1] )
			{
				//parallel
				case 'p':  
					args->parallel = 1;
					break;

				//k
				case 'k':  
					k = atoi(argv[i+1]);
					if(k == 0) {printf("k has to be non zero \n"); exit(-1); }
					args->k = k ;
					break;
	
				//number of trials
				case 't':  
					t = atoi(argv[i+1]);
					if(t == 0) {printf("t has to be non zero \n"); exit(-1); }
					args->trials = t ;
					break;
				
				case 'w':  
					w = atoi(argv[i+1]);
					if(w < 1 || w > 2) {printf("w has to be 1 or 2 \n"); exit(-1); }
					args->winSize = w ;
					break;
				
				case 'g':  
					
					if(argv[i][2] == 'b')
					{
						gbool = atoi(argv[i+1]);
						args->fillGaps = gbool;
					}
					else if(argv[i][2] == 's')
					{
						gsize = atoi(argv[i+1]);
						args->gapSize = gsize;
					}
					else if(argv[i][2] == 'u')
					{
						gu = atoi(argv[i+1]);
						args->GU = gu;
					}
					break;

				case 'f':
					printf("here: %s \n", argv[i+1]);
					strcpy(args->fileName, argv[i+1]);
				
				case 'r':
					printf("rnaup out: %s \n", argv[i+1]);
					strcpy(args->rnaupOut, argv[i+1]);
				
				default: 
					;
			} 
	}

	printf("Using k = %d \n", args->k);
	printf("Using parallel = %d \n", args->parallel);
	printf("using trials  = %d \n", args->trials);
	printf("using gaps  = %d \n", args->fillGaps);

	return args;
}


int getInput(char ** strArray1)
{
	int numOfRNA = readAndTokenizeRNAs(strArray1);
	
	/** 2. Get the (even|odd) type of each RNA **/
	tokenizeNames(numOfRNA);
	
	
	readEvenOdd();

	return numOfRNA;
}