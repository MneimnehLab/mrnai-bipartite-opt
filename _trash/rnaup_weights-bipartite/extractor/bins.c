#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include "bins.h"
#include "rnaSeq.h"


extern OrigRNASeq * origRNASequences[10];
extern InterDims interDimArray[10];

void destroyStruct_(EOB * n, int i);

//Even_Odd_Barrier :- EOB

unsigned long long rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}


Pointer * createEmptyStructure()
{
	Pointer * head = (Pointer*) malloc(sizeof(Pointer));
	
	head->elem = NULL;
	
	return head;
}

Pointer * makeSomeBinStruct(OrigRNASeq ** origRNASequences, int numOfRNA)
{
	Pointer * head = createEmptyStructure();
	
	int i;
	
	EOB * n = (EOB*) malloc(sizeof(EOB));
	n->next = NULL;
	n->origId = origRNASequences[0]->origId;
	n->type = origRNASequences[0]->type;
		
	head->elem = n;
	
	//printf("id: %d \n", head->elem->origId);
	for(i=1;i<numOfRNA;i++)
	{
		EOB * n2 = (EOB*)malloc(sizeof(EOB));
		n2->next = NULL;
		n2->origId = origRNASequences[i]->origId;
		n2->type = origRNASequences[i]->type;
		
		//printf("id: %d \n", n2->origId);
		
		n->next = n2;
		
		n = n2;
	}
	
	return head;
}

Pointer * makeRandomBinStruct(OrigRNASeq ** origRNASequences, int numOfRNA)
{
	//srand(time(NULL));
	srand(rdtsc());
	int array[MAX_NUM_RNA];
	int i;
	
	for(i=0;i<numOfRNA;i++)
	{
		array[i] = origRNASequences[i]->origId;
	}
	//now randomize it
	for(i=0;i<numOfRNA;i++)
	{
		//select index 2
		int r1 = rand() % numOfRNA;	//0,1,2,..,n-1

		int temp = array[i];
		array[i] = array[r1];
		array[r1] = temp;
	}
	
	
	Pointer * head = createEmptyStructure();
	
	EOB * n = (EOB*) malloc(sizeof(EOB));
	n->next = NULL;
	n->origId = origRNASequences[array[0]]->origId;
	//printf("n->origId = %d \n", n->origId);
	n->type = origRNASequences[array[0]]->type;
		
	head->elem = n;
	
	//printf("id: %d \n", head->elem->origId);
	for(i=1;i<numOfRNA;i++)
	{
		EOB * n2 = (EOB*) malloc(sizeof(EOB));
		n2->next = NULL;
		n2->origId = origRNASequences[array[i]]->origId;
		//printf("n2->origId = %d \n", n2->origId);
		n2->type = origRNASequences[array[i]]->type;
		
		//printf("id: %d \n", n2->origId);
		
		n->next = n2;
		
		n = n2;
	}
	
	return head;
}



Pointer * makeStructFromNums(char * numbers)
{
	//assume there will be 8, and that 0...3 are even, and 4...7 are odd
	int i;
	
	Pointer * head = createEmptyStructure();
	EOB * array[8];
	EOB * n, * last = NULL; 
	
	for(i=0;i<8;i++)
	{
		char c = numbers[i];
		
		n = (EOB*) malloc(sizeof(EOB));
		int num = atoi(&c);

		n->origId = num;
//printf(":: num = %d \n", num);
		if(num <= 3)
			n->type = 0;
		else
			n->type = 1;

		if(i != 0)
			last->next = n;
		else
			head->elem = n;

		last = n;
		
	}
	last->next = NULL;

	return head;
}

// EOB * flip_(EOB * curr, EOB * prev, int num, int stop)
// {
// 	if(num < stop)
// 		flip_(curr->next, curr, num + 1,  stop);
	
// 	curr->next = prev;
// }

void flipSubChain(Pointer * head, int start, int end)
{
	EOB * n = head->elem;
	EOB * startEOB, * endEOB;
	
	int i = 0;
	EOB * array[20];
	EOB * prevToStart = NULL;
	int flipCount = 0;
	while(n != NULL)
	{	
		if(i == start - 1)
			prevToStart = n;
			
		if(i == start)
		{
			startEOB = n;
			printf("Start id = %d \n", startEOB->origId);
		}	
		else if(i == end)
		{
			endEOB = n;
			printf("End id = %d \n", endEOB->origId);
		}
		
		if(i >= start && i <= end)
		{
			//printf("reversing pointers of %d \n", n->origId);
			
			array[flipCount++] = n;
		}
			
		n = n->next;
		i++;
	}
	
	EOB * nextToEnd = endEOB->next;
	
	for(i=flipCount-1;i>0;i--)
	{
		EOB * thisOne = array[i];
		EOB * realPrev = array[i-1];
		
		thisOne->next = realPrev;
	}
	
	startEOB->next = nextToEnd;
	if(prevToStart == NULL)	//, then prev is head pointer;
	{
		head->elem = endEOB;
	}
	else
		prevToStart->next = endEOB;
	
	
	printf("reversed... \n");
	//flip_(startEOB->next, startEOB, 0, end - start);
	
	//		else if(i > start && i < end)
//		{
			//flip_(n->next, n, num, int stop)
		//}
	/*
	//now that we have the sub chain, and have reversed the inside, just reverse end points
	EOB * prevOfStart, *prevOfEnd, * nextToEnd, *nextToStart;	
	
	prevOfStart = startEOB->prev;
	nextToStart = startEOB->next;
	
	prevOfEnd = endEOB->prev;
	nextToEnd = endEOB->next;
	
	startEOB->next = nextToEnd;
	endEOB->prev = prevOfStart;
	
	startEOB->prev = nextToStart;
	endEOB->next = prevOfEnd;
	
	if(nextToEnd != NULL)
	{
		printf("next to end not null \n");
		nextToEnd->prev = startEOB;
	}
	else
		printf("nextToEnd is null \n");
	
	if(prevOfStart != NULL) 
		prevOfStart->next = endEOB;
	else //prev was head pointer!
	{
		printf("here 1\n");
		head->elem = endEOB;
		if(endEOB == NULL)
			printf("endEOB is null \n");
		else
			printf("endEOB is not null \n");
			
			
		if(endEOB->prev != NULL)
		{	
			printf("endEOBV->prev is not null \n");
			printf("endEOB->prev = %d \n", endEOB->prev->origId);
		}
		else
			printf("endEOB->prev is null \n", endEOB->prev->origId);
		
		
		endEOB->prev = NULL;
	}
	if(nextToStart != NULL) nextToStart->prev = startEOB;
	
	//printChain(head);
	
	//if they are same type:
	// 		if they even or odd, ignore, else if barrier, get rid of it
	
	
	//printChain(head);
	//printf("\n\n");
	
	*/
}


/*
void flipSubChain(Pointer * head, int start, int end)
{
	EOB * n = head->elem;
	EOB * startEOB, * endEOB;
	
	int i = 0;
	while(n != NULL)
	{	
		EOB * hypoNext = n->next;
		int nextId =  (hypoNext != NULL) ? hypoNext->origId : -1;
		//printf("on %d, setting next to %d \n", n->origId, nextId);
		
		if(i == start)
		{
			startEOB = n;
			//printf("Start id = %d \n", startEOB->origId);
		}	
		else if(i == end)
		{
			endEOB = n;
			//printf("End id = %d \n", endEOB->origId);
		}
		//exchange next/prev of each item within {(_,_), not [_,_] } the list
		else if(i > start && i < end)
		{
			//printf("reversing pointers of %d \n", n->origId);
			EOB * temp = n->next;
			n->next = n->prev;
			n->prev = temp;
		}
			
		n = hypoNext;
		i++;
	}
	
	//now that we have the sub chain, and have reversed the inside, just reverse end points
	EOB * prevOfStart, *prevOfEnd, * nextToEnd, *nextToStart;	
	
	prevOfStart = startEOB->prev;
	nextToStart = startEOB->next;
	prevOfEnd = endEOB->prev;
	nextToEnd = endEOB->next;
	
	startEOB->next = nextToEnd;
	endEOB->prev = prevOfStart;
	
	startEOB->prev = nextToStart;
	endEOB->next = prevOfEnd;
	
	prevOfStart->next = endEOB;
	nextToStart->prev = startEOB;
	
	//Now see whats going on at end points
	//Only need to care about New Starts' prev and New End's next
	EOB * newStart = endEOB;
	EOB * newEnd = startEOB;
	
	printf("New Start type = %d \n", newStart->type);
	printf("New Start prev type = %d \n", newStart->prev->type);
	
	printf("New End type = %d \n", newEnd->type);
	printf("New End next type = %d \n", newEnd->next->type);
	
	printChain(head);
	
	//if they are same type:
	// 		if they even or odd, ignore, else if barrier, get rid of it
	if(newStart->type == 2 && newStart->prev->type == 2)
	{	
		//now we have to see what is on either side of the two barriers;
		//EOB * nextToStart = newStart->next;
		EOB * firstBarrier = newStart->prev;
		EOB * secondBarrier = newStart;
		
		printf("%d || %d \n", firstBarrier->prev->origId, secondBarrier->next->origId);
		
		if(firstBarrier->prev->type == secondBarrier->next->type)
		{
			//same type, delete both barriers
			printf("same type \n");
			firstBarrier->prev->next = secondBarrier->next;
			secondBarrier->next->prev = firstBarrier->prev;
			
			free(firstBarrier);
			free(secondBarrier);
		}
		else
		{
			//delete one barrier. choose start->prev, ie. firstBarrier
			printf("different type \n");
			firstBarrier->prev->next = secondBarrier;
			secondBarrier->prev = firstBarrier->prev;
		}
		
	}
	
	printChain(head);
	
	if(newEnd->type == 2 && newEnd->next->type == 2)
	{	
		//now we have to see what is on either side of the two barriers;
		EOB * firstBarrier = newEnd;
		EOB * secondBarrier = newEnd->next;
		
		printf("%d || %d \n", firstBarrier->prev->origId, secondBarrier->next->origId);
		
		if(firstBarrier->prev->type == secondBarrier->next->type)
		{
			//same type, delete both barriers
			printf("same type \n");
			firstBarrier->prev->next = secondBarrier->next;
			secondBarrier->next->prev = firstBarrier->prev;
			
			free(firstBarrier);
			free(secondBarrier);
		}
		else
		{
			//delete one barrier. choose start->prev, ie. firstBarrier
			printf("different type \n");
			firstBarrier->prev->next = secondBarrier;
			secondBarrier->prev = firstBarrier->prev;
		}
	}
	
	
	printChain(head);
	
	//if they are different types, introduce barrier
	if(newStart->type != newStart->prev->type)
	{
		EOB * barrier = malloc(sizeof(EOB));
		EOB * prevToStart = newStart->prev;
		
		barrier->next = newStart;
		barrier->prev = prevToStart;
		barrier->type = 2;
		barrier->origId = 987;
		
		newStart->prev = barrier;
		prevToStart->next = barrier;
	}
	if(newEnd->type != newEnd->next->type)
	{
		EOB * barrier = malloc(sizeof(EOB));
		EOB * nextToEnd = newEnd->next;
		
		barrier->next = nextToEnd;
		barrier->prev = newEnd;
		barrier->type = 2;
		barrier->origId = 345;
		
		newEnd->next = barrier;
		nextToEnd->prev = barrier;
	}
	
	/*
	
	
	printf("newStart->type = %d \n", newStart->type);
	printf("newStart->prev->type = %d \n", newStart->prev->type);
	//if they are different types, introduce barrier
	if(newStart->type != newStart->prev->type)
	{
		EOB * barrier = malloc(sizeof(EOB));
		EOB * prevToStart = newStart->prev;
		
		barrier->next = newStart;
		barrier->prev = prevToStart;
		barrier->type == 2;
		
		newStart->prev = barrier;
		prevToStart->next = barrier;
	}
	if(newEnd->type != newEnd->next->type)
	{
		EOB * barrier = malloc(sizeof(EOB));
		EOB * nextToEnd = newEnd->next;
		
		barrier->next = nextToEnd;
		barrier->prev = newEnd;
		barrier->type == 2;
		
		newEnd->next = barrier;
		nextToEnd->prev = barrier;
	}
	* /
	printChain(head);
	printf("\n\n");
}
*/


void printChain(Pointer * head)
{

	EOB * n = head->elem;
	int c = 0;
	printf("Chain: \n");
	while(n != NULL)
	{	
		if(n->type == 2)
			printf("| ");
		else if(n->type == 0)
			printf("e%d ", n->origId);
		else if(n->type == 1)
			printf("o%d ", n->origId);
			
		n = n->next;
		c++;
	}
	printf("\n\n");
}


/*
int main()
{
	Pointer * head = makeFirstBinStruct();
	printChain(head);

//	flipSubChain(head, 14, 26);
//	flipSubChain(head, 15, 26);
	flipSubChain(head, 13, 26);
	
	return 0;
}
*/


void determineStruct(Pointer * head)
{
	EOB * n = head->elem;
	
	int currKnown = n->type;
	int level = 1;
	printf("Chain: \n");
	printf("Level 1: ");
	while(n != NULL)
	{	
		currKnown = n->type;
	
		if(n->type == 2)
			printf("| ");
		else if(n->type == 0)
			printf("e%d ", n->origId);
		else if(n->type == 1)
			printf("o%d ", n->origId);
			
		n = n->next;
		
		if(n != NULL)
			if(currKnown != n->type)
			{
				printf("\nLevel %d: ", ++level);
				//printf("\n");
			}
	}
	printf("\n\n");

}

int getBBHeight(Pointer * head)
{
	EOB * n = head->elem;
	int currKnown = n->type;
	int height = 1;
	while(n != NULL)
	{	
		currKnown = n->type;
		n = n->next;
		
		if(n != NULL)
			if(currKnown != n->type)
			{
				height++;
			}
	}
	return height;
}
/*
void convertToTree(Pointer * head)
{
	Node * root = (Node *) createRoot(-1);
	
	EOB * n = head->elem;
	int currKnown = n->type;
	Node * node;
	Node * parent = root;
	while(n != NULL)
	{	
		currKnown = n->type;
		
		node = (Node *) createNode (n->origId, root);
		
		n = n->next;
		
		if(n != NULL)
			if(currKnown != n->type)
			{
				parent = n;
			}
	}
	
	
}
*/

SingleRunConfig * generateSubsetConfigFromBins(Pointer * head, int from, int to)
{
	EOB * n0 = head->elem;

	//printf("n0's origId: %d \n", n0->origId);		//printf("head->elem's origId: %d \n", head->elem->origId);
	int currKnown = n0->type;
	int level = 1;
	EOB * fromNode = NULL;
	while(n0 != NULL)
	{	

		currKnown = n0->type;

		if(level == from)
		{
			fromNode = n0;
			break;
		}

		n0 = n0->next;

		if(n0 != NULL)
			if(currKnown != n0->type)
				++level;

	}
	
	SingleRunConfig * data = (SingleRunConfig*)malloc(sizeof(SingleRunConfig));

	int i,j;
	//int setup[10][10];
	//int parentOf[10][10];
	int newId2OldId[10];
	int oldId2NewId[10];
	char ** strArray1;

	for(i=0;i<10;i++) data->newId2OldId[i] = -9;
	for(j=0;j<10;j++){newId2OldId[j] = -9;oldId2NewId[j] = -9;}
	for(i=0;i<10;i++) for(j=0;j<10;j++)	data->parentOf[i][j] = -9;
	for(i=0;i<10;i++) for(j=0;j<10;j++)	data->setup[i][j] = -9;
	
	//to create setup[][] array, just go thorugh the sequences linearly
	int num = 0;
	//for(i=from;i<=to;i++)
	EOB * n = n0;
	
	
	int setupRow = 0;
	int setupColumn = 0;
	
	char stringRep[100] = "";
	
	//printf("PP stringRep = %s \n", stringRep);
	while(level <= to && n != NULL)
	{
		currKnown = n->type;
		
		newId2OldId[num] = n->origId;
		data->newId2OldId[num] = n->origId;

		oldId2NewId[n->origId] = num;

		char t[10];
		sprintf(t,"%d",n->origId);	// converts to decimal base.
		
		strcat(stringRep, t);
		
		data->setup[setupRow][setupColumn] = num;

		data->rnaSequences[num] = (RNASeq*)malloc(sizeof(RNASeq));
		data->rnaSequences[num]->origId = n->origId;

		data->rnaSequences[num]->id = num;

		data->rnaSequences[num]->type = origRNASequences[n->origId]->type;
		data->rnaSequences[num]->string = origRNASequences[n->origId]->string;
		data->rnaSequences[num]->originalLength = origRNASequences[n->origId]->originalLength;

		strcpy(data->rnaSequences[num]->name, origRNASequences[n->origId]->name);

		num++;
		n = n->next;

		setupColumn++;
		if(n != NULL)
		{
			if(currKnown != n->type)
			{
				++level;
				setupRow++;
				setupColumn = 0;
				
				if(level <= to)
					strcat(stringRep, "|");
			}
			else
				strcat(stringRep, "_");
		}
	}
	/*
	//to get parent relations, we will have to traverse it in psuedo tree-style
	//also, this should be done only after newId2OldId & oldId2NewId maps are done both ways completely
	for(i=from;i<to;i++)                   
	{                       
           
		RowHead * head = heads[i];
		Node * n = head->first;
		
		while(n != NULL)
		{
			//look at n's children
			for(j=0;j<n->childCount;j++)
			{
				data->parentOf[oldId2NewId[n->id]][oldId2NewId[n->children[j]->id]] = 1;
			}
			
			n = n->next;
		}

	}
	*/

	data->numOfRNA = num;
	data->totalLevels = from - to + 1;
	//printf("data->totalLevels = %d \n", data->totalLevels);

	data->totalLevels = to - from + 1;
	//printf("data->totalLevels = %d \n", data->totalLevels);
	
	data->stringRep = (char*)malloc(sizeof(char) * (strlen(stringRep) + 1));
	strcpy(data->stringRep, stringRep);
	
	//printf("data->stringRep = %s \n", data->stringRep);
	//printf("1 stringRep = %s \n", stringRep);
	
	//fix parentOf matrix
	for(i=0;i<9;i++)	//rows
	{
		for(j=0;j<10;j++)	//cols on top level
		{
			if(data->setup[i][j] == -9) continue;
			
			int parentRNA = data->setup[i][j];
			int k;
			for(k=0;k<10;k++)	//cols on second level
			{
				if(data->setup[i+1][k] == -9) continue;
				
				int childRNA = data->setup[i+1][k];
				data->parentOf[parentRNA][childRNA] = 1;
				
				printf("setting parent=1 for %d > %d \n", parentRNA, childRNA);
			}
		}
	}
	printf("exiting...\n");

	
	return data;
}

void printFlatStruct2File(Pointer * head, FILE * file)
{
	EOB * n = head->elem;
	
	int currKnown = n->type;
	int level = 1;
	while(n != NULL)
	{	
		currKnown = n->type;
	
		if(n->type == 0)
			fprintf(file, "e%d ", n->origId);
		else if(n->type == 1)
			fprintf(file, "o%d ", n->origId);
			
		n = n->next;
		
		if(n != NULL)
			if(currKnown != n->type)
			{
				fprintf(file, "| ");
			}
	}
}

void printFlatStruct(Pointer * head)
{
	EOB * n = head->elem;
	
	int currKnown = n->type;
	int level = 1;
	while(n != NULL)
	{	
		currKnown = n->type;
	
		if(n->type == 0)
			printf("e%d ", n->origId);
		else if(n->type == 1)
			printf("o%d ", n->origId);
			
		n = n->next;
		
		if(n != NULL)
			if(currKnown != n->type)
			{
				printf("| ");
			}
	}
}

void destroyStruct(Pointer * head)
{
	EOB * n = head->elem;
	
	destroyStruct_(n, 0);
}


void destroyStruct_(EOB * n, int i)
{
	if(n->next == NULL)
	{
		free(n);
		return;
	}
	
	destroyStruct_(n->next, i+1);
	free(n);
}




void printReverseChain(Pointer * head)
{
/*
	EOB * n = head->elem;
	
	while(1)
	{
		if(n->next == NULL)
			break;
			
		n = n->next;
	}
	printf("Reverse Chain: \n");
	
	while(n != NULL)
	{	
		if(n->type == 2)
			printf("| ");
		else if(n->type == 0)
			printf("e%d ", n->origId);
		else if(n->type == 1)
			printf("o%d ", n->origId);
			
		n = n->prev;
	}
	printf("\n\n");
	*/
}










































































/*
Pointer * makeFirstBinStruct()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n = (EOB*) malloc(sizeof(EOB));
	n->next = NULL;
	n->origId = 0;
	n->type = 0;
	
	EOB * n3 = (EOB*) malloc(sizeof(EOB));
	n3->next = NULL;
	n3->type = 2;
	n3->origId = 99;
	
	n->next = n3;
	
	head->elem = n;
	n = n3;
	
	int next = 1;
	for(i=1;i<20;i++)
	{
		EOB * n2;
		if(next == 0)
		{
			n2 = (EOB*) malloc(sizeof(EOB));
			n2->next = NULL;
			n2->origId = i;
			n2->type = 0;
	
			next = 1;
		}
		else
		{
			n2 = (EOB*) malloc(sizeof(EOB));
			n2->next = NULL;
			n2->origId = i;
			n2->type = 1;
	
			next = 0;
		}
		
		n->next = n2;
		
		if(i != 20 - 1)
		{
			EOB * n3 = (EOB*) malloc(sizeof(EOB));
			n3->next = NULL;
			n3->type = 2;
			n3->origId = 100 + i;
			
			n2->next = n3;
			
			n = n3;
		}
		
		printf("%d ", n2->type);
	}
	printf("\n");
	
	
	return head;
}
*/


/*
Pointer * makeSampleBinStruct1()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = (EOB*) malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = (EOB*) malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = (EOB*) malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = (EOB*) malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = (EOB*) malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = (EOB*) malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = (EOB*) malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = (EOB*) malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n1->next = n2;
	n2->next = n3;
	n3->next = n7;
	n7->next = n6;
	n6->next = n0;
	n0->next = n4;
	n4->next = n5;
	n5->next = NULL;
	
	head->elem = n1;
	
	
	
	return head;
}

Pointer * makeSampleBinStruct2()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = (EOB*) malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = (EOB*) malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = (EOB*) malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = (EOB*) malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = (EOB*) malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = (EOB*) malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = (EOB*) malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = (EOB*) malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n6->next = n4;
	n4->next = n7;
	n7->next = n5;
	n5->next = n3;
	n3->next = n2;
	n2->next = n0;
	n0->next = n1;
	n1->next = NULL;
	
	head->elem = n6;
	
	
	
	return head;
}

Pointer * makeSampleBinStruct3()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = (EOB*) malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = (EOB*) malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = (EOB*) malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = (EOB*) malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n4->next = n6;
	n6->next = n5;
	n5->next = n7;
	n7->next = n2;
	n2->next = n0;
	n0->next = n3;
	n3->next = n1;
	n1->next = NULL;
	
	head->elem = n4;
	
	
	
	return head;
}

Pointer * makeSampleBinStruct4()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n1->next = n6;
	n6->next = n0;
	n0->next = n2;
	n2->next = n5;
	n5->next = n7;
	n7->next = n4;
	n4->next = n3;
	n3->next = NULL;
	
	head->elem = n1;
	
	
	
	return head;
}


Pointer * makeSampleBinStruct5()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n2->next = n4;
	n4->next = n7;
	n7->next = n6;
	n6->next = n3;
	n3->next = n0;
	n0->next = n5;
	n5->next = n1;
	n1->next = NULL;
	
	head->elem = n2;
	
	
	
	return head;
}


Pointer * makeSampleBinStruct6()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n0->next = n5;
	n5->next = n6;
	n6->next = n7;
	n7->next = n1;
	n1->next = n3;
	n3->next = n2;
	n2->next = n4;
	n4->next = NULL;
	
	head->elem = n0;
	
	
	
	return head;
}



Pointer * makeSampleBinStruct15()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n6->next = n0;
	n0->next = n2;
	n2->next = n5;
	n5->next = n7;
	n7->next = n4;
	n4->next = n1;
	n1->next = n3;
	n3->next = NULL;
	
	head->elem = n6;
	
	
	
	return head;
}


Pointer * makeSampleBinStruct16()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n7->next = n2;
	n2->next = n5;
	n5->next = n1;
	n1->next = n3;
	n3->next = n4;
	n4->next = n0;
	n0->next = n6;
	n6->next = NULL;
	
	head->elem = n7;
	
	
	
	return head;
}


Pointer * makeSampleBinStructCopACopT()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 1;
	
	n0->next = n1;
	n1->next = NULL;
	
	head->elem = n0;
	
	
	
	return head;
}


Pointer * makeSampleBinStructJustTwo()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 1;
	
	n0->next = n1;
	n1->next = NULL;
	
	head->elem = n0;
	
	
	
	return head;
}

Pointer * makeU6_U4U2()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 1;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 1;
	
	n0->next = n1;
	n1->next = n2;
	n2->next = NULL;
	
	head->elem = n0;
	
	
	
	return head;
}

Pointer * makeSampleBinStructJustTwoStartOdd()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 1;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	n0->next = n1;
	n1->next = NULL;
	
	head->elem = n0;
	
	
	
	return head;
}

Pointer * justFourOddToEven()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 1;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 1;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	n0->next = n1;
	n1->next = n2;
	n2->next = n3;
	n3->next = NULL;
	
	head->elem = n0;
	
	
	
	return head;
}



Pointer * justFourOddToEven2_31_0()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 1;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 1;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	n2->next = n3;
	n3->next = n1;
	n1->next = n0;
	n0->next = NULL;
	
	head->elem = n2;
	
	
	
	return head;
}



Pointer * justFourOddToEven3_2_1_0()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 1;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 1;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	n3->next = n2;
	n2->next = n1;
	n1->next = n0;
	n0->next = NULL;
	
	head->elem = n3;
	
	
	
	return head;
}


Pointer * makeSampleBinStruct17()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;

	//Level 1: e2 e3
	//Level 2: o4 o5
	//Level 3: e1
	//Level 4: o7 o6
	//Level 5: e0
	
	
	n2->next = n3;
	n3->next = n4;
	n4->next = n5;
	n5->next = n1;
	n1->next = n7;
	n7->next = n6;
	n6->next = n0;
	n0->next = NULL;
	
	head->elem = n2;
	
	
	
	return head;
}


Pointer * testConvgOutStruct1()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	//o5 | e1 e3 e2 | o4 o7 | e0 | o6
	
	n5->next = n1;
	n1->next = n3;
	n3->next = n2;
	n2->next = n4;
	n4->next = n7;
	n7->next = n0;
	n0->next = n6;
	n6->next = NULL;
	
	head->elem = n5;
	
	
	
	return head;
}










Pointer * somethingWrongInCnvg1()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n1->next = n7;
	n7->next = n0;
	n0->next = n4;
	n4->next = n6;
	n6->next = n3;
	n3->next = n2;
	n2->next = n5;
	n5->next = NULL;
	
	head->elem = n1;
	
	
	
	return head;
}


Pointer * somethingWrongInCnvg2()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 1;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 1;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 0;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 0;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n5->next = n2;
	n2->next = n3;
	n3->next = n6;
	n6->next = n7;
	n7->next = n1;
	n1->next = n4;
	n4->next = n0;
	n0->next = NULL;
	
	head->elem = n5;
	
	
	
	return head;
}



Pointer * testMixedConvg1()
{
	//o1 o7 o5 | e0 e6 e2 | o3 | e4
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 1;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 1;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 0;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 0;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n1->next = n7;
	n7->next = n5;
	n5->next = n0;
	n0->next = n6;
	n6->next = n2;
	n2->next = n3;
	n3->next = n4;
	n4->next = NULL;
	
	head->elem = n1;
	
	
	
	return head;
}


Pointer * createFromString(char * str)
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * node, * node2;
	printf("%c \n", str[0]);
	int n = atoi(str[0]);
	printf("d:  %d \n", n);
	node = malloc(sizeof(EOB));
	printf("A\n");
	node->origId = n;
	printf("B\n");
	node->type = (n % 2 == 0) ? 0 : 1;
		printf("C\n");
	head->elem = node;
	printf("D\n");
	for(i=1;i<8;i++)
	{
		printf("%c \n", str[i]);
		n = atoi(str[i]);
		printf("d:  %d \n", n);
		node2 = malloc(sizeof(EOB));
		node2->origId = n;
		node2->type = (n % 2 == 0) ? 0 : 1;
		
		node->next = node2;
		
		node = node2;
	}
	
	node->next = NULL;
	
	return head;
}














Pointer * somethingWrongInCnvg_16()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n2->next = n0;
n0->next = n7;
n7->next = n6;
n6->next = n3;
n3->next = n1;
n1->next = n4;
n4->next = n5;
n5->next = NULL;
head->elem = n2;
	
	
	
	return head;
}


Pointer * somethingWrongInCnvg_14()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
	n6->next = n0;
n0->next = n7;
n7->next = n4;
n4->next = n5;
n5->next = n2;
n2->next = n1;
n1->next = n3;
n3->next = NULL;
head->elem = n6;
	
	
	return head;
}

Pointer * somethingWrongInCnvg_12()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
n4->next = n5;
n5->next = n7;
n7->next = n6;
n6->next = n3;
n3->next = n1;
n1->next = n2;
n2->next = n0;
n0->next = NULL;
head->elem = n4;

	
	return head;
}



Pointer * somethingWrongInCnvg_17()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
n1->next = n3;
n3->next = n5;
n5->next = n0;
n0->next = n6;
n6->next = n2;
n2->next = n4;
n4->next = n7;
n7->next = NULL;
head->elem = n1;
	
	return head;
}





Pointer * somethingWrongInMixedCnvg_1()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 1;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 1;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 0;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 0;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
n6->next = n4;
n4->next = n0;
n0->next = n5;
n5->next = n7;
n7->next = n3;
n3->next = n1;
n1->next = n2;
n2->next = NULL;
head->elem = n6;
	
	return head;
}


Pointer * mix11()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 1;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 1;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 0;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 0;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	//o7 | e4 e6 e2 e0 | o5 o3 o1

n7->next = n4;
n4->next = n6;
n6->next = n2;
n2->next = n0;
n0->next = n5;
n5->next = n3;
n3->next = n1;
n1->next = NULL;
head->elem = n7;
	
	return head;
}



Pointer * somethingWrongInMixedCnvg_2()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
n1->next = n7;
n7->next = n5;
n5->next = n0;
n0->next = n6;
n6->next = n2;
n2->next = n3;
n3->next = n4;
n4->next = NULL;
head->elem = n1;
	
	return head;
}



Pointer * somethingWrongInCnvg_20()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
n2->next = n1;
n1->next = n3;
n3->next = n7;
n7->next = n4;
n4->next = n5;
n5->next = n0;
n0->next = n6;
n6->next = NULL;
head->elem = n2;
	
	return head;
}



Pointer * somethingWrongInCnvg_B_15()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	
n3->next = n5;
n5->next = n0;
n0->next = n6;
n6->next = n2;
n2->next = n7;
n7->next = n1;
n1->next = n4;
n4->next = NULL;
head->elem = n3;
	
	return head;
}


Pointer * somethingWrongInCnvg_B_4()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
		
	n0->next = n1;
	n1->next = n6;
	n6->next = n4;
	n4->next = n2;
	n2->next = n7;
	n7->next = n3;
	n3->next = n5;
	n5->next = NULL;
	head->elem = n0;
		
	return head;
}


Pointer * somethingWrongInCnvg_B_10()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
		
	n6->next = n0;
	n0->next = n1;
	n1->next = n4;
	n4->next = n5;
	n5->next = n7;
	n7->next = n3;
	n3->next = n2;
	n2->next = NULL;
	head->elem = n6;
		
	return head;
}










Pointer * makeBinStruct_Highest()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	//e0 | o6 | e2 e1 | o7 o5 | e3 | o4 
	
	n0->next = n6;
	n6->next = n2;
	n2->next = n1;
	n1->next = n7;
	n7->next = n5;
	n5->next = n3;
	n3->next = n4;
	n4->next = NULL;
	
	head->elem = n0;
	
	return head;
}

Pointer * makeBinStruct_2ndHighest()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
	//e0 | o6 | e3 | o5 o4 | e2 e1 | o7 

	
	n0->next = n6;
	n6->next = n3;
	n3->next = n5;
	n5->next = n4;
	n4->next = n2;
	n2->next = n1;
	n1->next = n7;
	n7->next = NULL;
	
	head->elem = n0;
	
	return head;
}



Pointer * runB_no6()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 0;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	EOB * n6 = malloc(sizeof(EOB));
	n6->origId = 6;
	n6->type = 1;
	
	EOB * n7 = malloc(sizeof(EOB));
	n7->origId = 7;
	n7->type = 1;
	
//e3 e1 | o4 o5 | e2 e0 | o7 o6

	
	n3->next = n1;
	n1->next = n4;
	n4->next = n5;
	n5->next = n2;
	n2->next = n0;
	n0->next = n7;
	n7->next = n6;
	n6->next = NULL;
	
	head->elem = n3;
	
	return head;
}





Pointer * nice6_1()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 1;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	
//o3 o4 | e1 e0 | o5 | e2  :: 

	
	n3->next = n4;
	n4->next = n1;
	n1->next = n0;
	n0->next = n5;
	n5->next = n2;
	n2->next = NULL;
	
	head->elem = n3;
	
	return head;
}

Pointer * nice6_2()
{
	int i;
	
	Pointer * head = createEmptyStructure();
	
	EOB * n0 = malloc(sizeof(EOB));
	n0->origId = 0;
	n0->type = 0;
	
	EOB * n1 = malloc(sizeof(EOB));
	n1->origId = 1;
	n1->type = 0;
	
	EOB * n2 = malloc(sizeof(EOB));
	n2->origId = 2;
	n2->type = 0;
	
	EOB * n3 = malloc(sizeof(EOB));
	n3->origId = 3;
	n3->type = 1;
	
	EOB * n4 = malloc(sizeof(EOB));
	n4->origId = 4;
	n4->type = 1;
	
	EOB * n5 = malloc(sizeof(EOB));
	n5->origId = 5;
	n5->type = 1;
	
	
//o5 | e2 | o4 | e0 | o3 | e1 

	
	n5->next = n2;
	n2->next = n4;
	n4->next = n0;
	n0->next = n3;
	n3->next = n1;
	n1->next = NULL;
	
	head->elem = n5;
	
	return head;
}
*/