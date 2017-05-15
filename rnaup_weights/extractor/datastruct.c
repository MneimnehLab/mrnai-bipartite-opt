#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <float.h> 
#include <stdarg.h>
#include "rnaseq.h"
#include "datastruct.h"

extern OrigRNASeq * origRNASequences[10];
extern InterDims interDimArray[10];

void getHeight_(Node * n, int level, int * maxH);

RowHead * heads[10];
int rowCount = 0;

struct Node
{
	int id; 
	int childCount;
	int isRoot;
	
	Node * next;
	Node * prev;
	Node * parent;
	Node * children[10];
	
	RowHead * rowHead;
	
	OrigRNASeq origSeq;
};

struct RowHead
{
	int id;
	Node * first;
};

void set1()
{
	//create a test scenario
	Node n1, n2, n3, n4, n5, n6;
	
	n1.id = 1;
	n2.id = 2;
	n3.id = 3;
	n4.id = 4;
	n5.id = 5;
	n6.id = 6;
	
	n1.parent = NULL;
	n2.parent = &n1;
	n3.parent = &n1;
	n4.parent = &n2;
	n5.parent = &n3;
	n6.parent = &n2;
	
}

void set2()
{
	Node * root = (Node *) createRoot(0);
	//For now build the tree Left recursively. Have to fix this later.
	//That is, middle insertion.
	//when printing, see the difference betweem 7,10,14,11 and 7,10,11,14
	/*
	Node * n1 = (Node *) createNode(1, root);
	Node * n2 = (Node *) createNode(2, root);
	Node * n3 = (Node *) createNode(3, n1);
	Node * n4 = (Node *) createNode(4, n2);
	Node * n5 = (Node *) createNode(5, n3);
	Node * n6 = (Node *) createNode(6, n2);
	Node * n7 = (Node *) createNode(7, n5);
	Node * n8 = (Node *) createNode(8, n6);
	Node * n9 = (Node *) createNode(9, n6);
	Node * n10 = (Node *) createNode(10, n5);
	Node * n11 = (Node *) createNode(11, n9);
	Node * n12 = (Node *) createNode(12, n7);
	Node * n13 = (Node *) createNode(13, n12);
	Node * n14 = (Node *) createNode(14, n5);
	*/
	
	Node * n1  = (Node *) createNode (1, root);
	Node * n3  = (Node *) createNode (3, n1  );
	Node * n5  = (Node *) createNode (5, n3  );
	Node * n7  = (Node *) createNode (7, n5  );
	Node * n12 = (Node *) createNode (12, n7 );
	Node * n13 = (Node *) createNode (13, n12);
	Node * n10 = (Node *) createNode (10, n5 );
	Node * n14 = (Node *) createNode (14, n5 );
	Node * n2  = (Node *) createNode (2, root);
	Node * n4  = (Node *) createNode (4, n2  );
	Node * n6  = (Node *) createNode (6, n2  );
	Node * n8  = (Node *) createNode (8, n6  );
	Node * n9  = (Node *) createNode (9, n6  );
	Node * n11 = (Node *) createNode (11, n9 );
	
	

	/*printf("n6->id = %d \n", n6->id);
	printf("n2->id = %d \n", n2->id);
	printf("n2->children[0]->id = %d \n", n2->children[1]->id);*/

	/*
	printTree(root);
	printf("total rows = %d \n", rowCount);
	int i;
	for(i=0;i<rowCount;i++)
	{
		printRow(i);
	}
	
	printf("getting location of n6 \n");
	getLocationAsChild(n6);
	*/
	
	printSubset(2,4);
	
//	generateSubsetConfig(2, 4);
}


Node * createRoot(int id)
{
	//first clear RowHead array
	rowCount = 0;
	int i;
	for(i=0;i<11;i++)
		heads[i] = NULL;


	//printf("creating root \n");
	Node * n = (Node*)malloc(sizeof(Node));
	n->id = id;
	n->childCount = 0;
	n->next = NULL;
	n->prev = NULL;
	n->isRoot = 1;
	
	//create row header
	//printf("creating root row Header\n");
	RowHead * head = (RowHead*) malloc(sizeof(RowHead));
	head->id = 0;
	head->first = n;
	n->rowHead = head;
	
	//printf("set heads[0] = head\n");
	heads[rowCount++] = head;
		
	//printf("created root \n");
	return n;
}

Node * createNode(int id, Node * parent)
{
	Node * n = (Node*)malloc(sizeof(Node));
	n->id = id;
	n->childCount = 0;
	n->parent = parent;
	n->isRoot = 0;
	
	parent->children[parent->childCount] = n;
	parent->childCount++;
	
	int parentHeadId = parent->rowHead->id;

	//if this guy is the first node on this row, then curr id of rowCount should be parent->head->id + 1
	if(rowCount == parentHeadId + 1)
	{
		//create head
		RowHead * head = (RowHead*)malloc(sizeof(RowHead));
		head->id = rowCount;
		n->rowHead = head;
		head->first = n;
		heads[rowCount++] = head;
		n->next = NULL;
		n->prev = NULL;
	}
	else	//head for this level already exists
	{
		n->rowHead = heads[parentHeadId + 1];

		(getLastInRow(n->rowHead->id))->next = n;
	}
	
	return n;
}

void insertInRow(Node * n)
{
	//assume we always put at the end of siblings
	//so insertion in middle happens when 
	
	Node * lastSibling = (Node *) getLastChild(n->parent);
	//if(lastSibling == NULL)
}

Node * getLastChild(Node * parent)
{
	if(parent->childCount == 0)
		return NULL;
		
	return parent->children[parent->childCount-1];
}

void printTree(Node * n)
{
	printf("%d \n", n->id);
	
	if(n->childCount == 0)
		return;
		
	int i;
	for(i=0;i<n->childCount;i++)
	{
		printTree(n->children[i]);
	}
}

void printFlatTree(Node * n)
{
	if(n->isRoot == 0)
		printf("%d", n->id);
	
	if(n->childCount == 0)
		return;
		
	int i;
	for(i=0;i<n->childCount;i++)
	{
		printf("(");
		printFlatTree(n->children[i]);
		printf(")");
	}
}

void printFlatTree2File(Node * n, FILE * file)
{
	if(n->isRoot == 0)
		fprintf(file, "%d", n->id);
	
	if(n->childCount == 0)
		return;
		
	int i;
	for(i=0;i<n->childCount;i++)
	{
		fprintf(file, "(");
		printFlatTree2File(n->children[i], file);
		fprintf(file, ")");
	}
}



void printRow(int rowId)
{
	
	Node * n = heads[rowId]->first;
	while(n != NULL)
	{
		printf("%d ", n->id);
		n = n->next;
	}
	printf("\n");
}

Node * getLastInRow(int rowId)
{
	Node * n = heads[rowId]->first;
	if(n == NULL) return NULL;
	while(n->next != NULL)
	{
		n = n->next;
	}
	return n;
}

int getHeight(Node * n)
{
	int maxH = 0;
	getHeight_(n, 0, &maxH);
	return maxH;
}

void getHeight_(Node * n, int level, int * maxH)
{
	if(n->childCount == 0)
	{
		if(level > *maxH)
			*maxH = level;
	}
		
	int i;
	for(i=0;i<n->childCount;i++)
	{
		getHeight_(n->children[i], level+1, maxH);
	}
}


void printSubset(int from, int to)
{
	int i,j;
	
	for(i=from;i<=to;i++)
	{
		RowHead * head = heads[i];
		Node * n = head->first;
		
		while(n != NULL)
		{
			printf("%d ", n->id);
			n = n->next;
		}
		printf("\n");
	}
}

int getLocationAsChild(Node * n)
{
	if(n->parent == NULL) return -1;
	
	if(n->parent->childCount == 0)
		return -1;
	
	
	int i;
	for(i=0;i<n->parent->childCount;i++)
		if(n == n->parent->children[i])
			return i;
	return -1;
}


SingleRunConfig * generateSubsetConfig(int from, int to)
{
	SingleRunConfig * data = (SingleRunConfig*) malloc(sizeof(SingleRunConfig));

	int i,j;
	//int setup[10][10];
	//int parentOf[10][10];
	int newId2OldId[10];
	int oldId2NewId[10];
	char ** strArray1;
	
	for(j=0;j<10;j++){newId2OldId[j] = -9;oldId2NewId[j] = -9;}
	for(i=0;i<10;i++) for(j=0;j<10;j++)	data->parentOf[i][j] = -9;
	for(i=0;i<10;i++) for(j=0;j<10;j++)	data->setup[i][j] = -9;
	
	//to create setup[][] array, just go thorugh the sequences linearly
	int num = 0;
	for(i=from;i<=to;i++)
	{

		RowHead * head = heads[i];
		Node * n = head->first;
		int k = 0;
		while(n != NULL)
		{
			//printf("Looking at n->id = %d \n", n->id);
			newId2OldId[num] = n->id;
			oldId2NewId[n->id] = num;
			
			data->setup[i-from][k++] = num;
			
			data->rnaSequences[num] = (RNASeq*)malloc(sizeof(RNASeq));
			data->rnaSequences[num]->origId = n->id;
			data->rnaSequences[num]->id = num;
			data->rnaSequences[num]->type = origRNASequences[n->id]->type;
			data->rnaSequences[num]->string = origRNASequences[n->id]->string;
			data->rnaSequences[num]->originalLength = origRNASequences[n->id]->originalLength;
			strcpy(data->rnaSequences[num]->name, origRNASequences[n->id]->name);
/*
			printf("Orig ID = %d \n", n->id);			
			printf("Num = %d \n", num);
			printf("Id = %d \n", data->rnaSequences[num]->id );
			printf("Type = %d \n", data->rnaSequences[num]->type );
			printf("String = %s \n", data->rnaSequences[num]->string );
			printf("Original Len = %d \n", data->rnaSequences[num]->originalLength );
*/			
			num++;
			n = n->next;
		}
		
	}
	
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
	
	data->numOfRNA = num;
	data->totalLevels = from - to + 1;
	//printf("data->totalLevels = %d \n", data->totalLevels);
	
	data->totalLevels = to - from + 1;
	//printf("data->totalLevels = %d \n", data->totalLevels);

	
	return data;
}

void deallocData(SingleRunConfig * data)
{
	int num = data->numOfRNA;
//	printf("deallocing data \n");
//	printf("num = %d \n", num);
	
//	printf("str1 = %s \n", data->rnaSequences[0]->string);
//	printf("str2 = %s \n", data->rnaSequences[1]->string);
	
	int i;
	for(i=0;i<num;i++)
	{
		free(data->rnaSequences[i]->compressedRNA);
//		printf("free %d - 2 \n", i);
		free(data->rnaSequences[i]->expandedRNAmap);
//		printf("free %d - 3 \n", i);
		free(data->rnaSequences[i]);
//		printf("free %d - 4 \n", i);
	}
	free(data);
//	printf("free 5 \n");
}


void swapNodePointer(Node ** p1, Node ** p2)
{
	Node ** temp = p1;
	p1 = p2;
	p2 = temp;
}

/** Written for linear trees only!! **/
void swapLevels(int l1, int l2)
{
//printf("rowCount = %d \n", rowCount);
int i;
//printf("Rows: ");
//for(i=0;i<rowCount;i++)
//	if(heads[i]->first != NULL) 
//		printf("%d ", heads[i]->first->id);
//	else
//		printf("X ");
//printf("\n");

//printf("SL 1 , l1 = %d, l2 = %d\n", l1, l2);
	Node * n1 = heads[l1]->first;
	Node * n2 = heads[l2]->first;
//if(n1 == NULL)
//	printf("n1 NULL \n");
//else
//	printf("n1->id = %d \n", n1->id);
//if(n2 == NULL)
//	printf("n2 NULL \n");
//else
//	printf("n2->id = %d \n", n2->id);

//printf("SL 2 \n");
	RowHead * tmpR = n1->rowHead;
//printf("SL 2 1\n");
	n1->rowHead = n2->rowHead;
//printf("SL 2 3\n");
	n2->rowHead = tmpR;
//printf("SL 3 \n");
	tmpR = heads[l1];
	heads[l1] = heads[l2];
	heads[l2] = tmpR;
//printf("SL 4 \n");	
	//swap childrens' parents of both
	Node * temp = n2->children[0]->parent;
	n2->children[0]->parent = n1->children[0]->parent;
	n1->children[0]->parent = temp;
//printf("SL 5 \n");
	
	//swap children of both
	temp = n2->children[0];
	n2->children[0] = n1->children[0];
	n1->children[0] = temp;
//printf("SL 6 \n");
	//swap parent's children of both
	temp = n2->parent->children[0];
	n2->parent->children[0] = n1->parent->children[0];
	n1->parent->children[0] = temp;
//printf("SL 7 \n");
	//swap parents of both
	temp = n2->parent;
	n2->parent = n1->parent;
	n1->parent = temp;
//printf("SL 8 \n");	
		
	/*
	Node * temp = n1->parent;
	n1->parent = n2->parent;
	n2->parent = temp;
	
	//change n2's child
	if(n2->childCount !=0)
	{
		temp = n2->children[0];
		n1->children[0] = n2->children[0];
		n2->children[0] = temp;
	}
	
	if(n2->childCount !=0)
	{
		temp = n2->children[0];
		n1->children[0] = n2->children[0];
		n2->children[0] = temp;
	}
	*/
	/*
	Node ** a = &(n1->parent);
	Node ** b = (n1->parent != NULL) ? &(n1->parent->children[0]) : NULL;
	Node ** c = (n1->childCount != 0) ? &(n1->children[0]->parent) : NULL;
	Node ** d = (n1->childCount != 0) ? &(n1->children[0]) : NULL;			
	Node ** i = &(n2->parent);
	Node ** j = (n2->parent != NULL) ? &(n2->parent->children[0]) : NULL;
	Node ** k = (n2->childCount != 0) ? &(n2->children[0]->parent) : NULL;
	Node ** l = (n2->childCount != 0) ? &(n2->children[0]) : NULL;		
	
	//RowHead * x = n1->head;
	//RowHead * y = n2->head;
	//n1->head;
	
	swapNodePointer(a,i);
	swapNodePointer(b,j);
	swapNodePointer(c,k);
	swapNodePointer(d,l);
	
	//RowHead * tempR = heads[l1];
	//heads[l1] = heads[l2];
	//heads[l2] = tempR;
	*/
}              



void writeNames(Node * n, FILE * file)
{
	if(n->id == -1)
		return;

	if(n->isRoot == 0)
	{
		printf("%s \n", origRNASequences[n->id]->name);
		fprintf(file, "%s ", origRNASequences[n->id]->name);
	}
	
	if(n->childCount == 0)
		return;
		
	int i;
	for(i=0;i<n->childCount;i++)
	{
		fprintf(file, "(");
		writeNames(n->children[i], file);
		fprintf(file, ")");
	}
}

void printRows()
{
	int i;
	printf("Rows: ");
	for(i=0;i<rowCount;i++)
		if(heads[i]->first != NULL) 
			printf("%d ", heads[i]->first->id);
		else
			printf("X ");
	printf("\n");
}

void printRows2File(FILE * file)
{
	int i;
	fprintf(file, "Rows: ");
	for(i=0;i<rowCount;i++)
		if(heads[i]->first != NULL) 
			fprintf(file, "%d ", heads[i]->first->id);
		else
			fprintf(file, "X ");
	fprintf(file, "\n");
}



//9 4 5 2 6 0 7 3 8 1
int breaker()
{
	if(	heads[1]->first->id == 4 &&
		heads[2]->first->id == 5 &&
		heads[3]->first->id == 2 &&
		heads[4]->first->id == 6 &&
		heads[5]->first->id == 0 &&
		heads[6]->first->id == 7 &&
		heads[7]->first->id == 3 &&
		heads[8]->first->id == 8 &&
		heads[9]->first->id == 1 ) return 1;
	else return 0;
}
