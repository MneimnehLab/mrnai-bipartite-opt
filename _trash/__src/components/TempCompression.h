#ifndef _TEMP_COMPRESS_H
#define _TEMP_COMPRESS_H


/**
In this module, we filter a nucleotide if it is "bad" for all windows ending on it.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <float.h> 
#include <stdarg.h>
#include "wins.h"
#include "rnaseq.h"

#define TEST 0
/**
Current we have rnaWins[r1_3'][r2_3'][w]. We need something to flip the indices of RNA because we will run from left to right.
Make it something like rnaWins[r1_3'][r2_5  '][w]
*/


class TempCompression
{

public:

TempCompression()
{
    compLen = 0;
}
SingleRunConfig * data2;

//int compLen = 0, compLen = 0;
int compLen;
int compressedRNAsTempHold[10][10][MAX_RNA_SIZE];   //compressedRNAsTempHold[x][y] where x = the compressed RNA - 1, y = found with - 1
int interactsWith[10][10];
int compressedRNAsLens[5];



// void doCompression(SingleRunConfig * data1, int num);
// void relations(int num);
// void merge(int num);
// void showGoodNucsOnEven(int even, int odd);
// void showGoodNucsOnOdd(int odd, int even);



void doCompression(SingleRunConfig * data1, int num)
{
    data2 = data1;
    //printf("string 1: %s \n", data2->rnaSequences[0]->string);
//  printf("Doing compression...() \n");

#if TEST == 1
    int x;
    for(x=0;x<num;x++)
        printf("string %d: %s \n", x, data2->rnaSequences[x]->name);
#endif

    /*
    printf("string 1: %s \n", data2->rnaSequences[0]->string);
    printf("string 2: %s \n", data2->rnaSequences[1]->string);
    
    printf("A string 1 len: %d \n", data2->rnaSequences[0]->originalLength);
    printf("A string 2 len: %d \n", data2->rnaSequences[1]->originalLength);
    */
    char ** str1;
    
    int i,j,w,k;
    //printf("here 1\n");
    for(i=0;i<10;i++)
        for(j=0;j<10;j++)
            interactsWith[i][j] = -9;
    //printf("here 2\n");       
    for(i=0;i<10;i++)
        for(j=0;j<10;j++)
            for(k=0;k<MAX_RNA_SIZE;k++)
                compressedRNAsTempHold[i][j][k] = -9;
    //printf("here 3\n");
    for(i=0;i<num;i++)
    {
        data2->rnaSequences[i]->compressedRNA = (int *) malloc(sizeof(int) * MAX_RNA_SIZE); 
        data2->rnaSequences[i]->expandedRNAmap = (int *) malloc(sizeof(int) * MAX_RNA_SIZE);
    }


    relations(num);

    merge(num);
    
//  printf("D string 1 len: %d \n", data2->rnaSequences[0]->originalLength);
//  printf("D string 2 len: %d \n", data2->rnaSequences[1]->originalLength);

#if TEST == 1
    for(i=0;i<num;i++)
    {
        printf("RNA: %d :-  ", i);
        printf("Comprssd / Original Len: %d / %d \n ", data2->rnaSequences[i]->compressedLength, data2->rnaSequences[i]->originalLength);
    }
#endif

}

void relations(int num)
{
#if TEST == 1
    printf("relations(), num = %d \n", num);
    
    int v;
    for(v=0;v<num;v++)
        printf("B string 1 : %s \n", data2->rnaSequences[v]->string);
#endif
    
    
    
    //for each RNA, find out which RNAs interact with it.
    //get locations of interactions
    //now determine if RNA is even or odd
    //if even, run showGoodNucsOnEven(current, partner)
    //else, run showGoodNucsOnOdd(current, partner)
    
    int rna, partner;
    
    for(rna=0;rna<num;rna++)
    {
//      printf("rna = %d \n", rna);
        for(partner=0;partner<num;partner++)
        {
//          printf("partner = %d \n", partner);
            if(!(data2->parentOf[rna][partner] == -9 && data2->parentOf[partner][rna] == -9))
            {
//              printf("A \n");
                if(data2->rnaSequences[rna]->type == 0)
                {
//                  printf("B \n");
                    showGoodNucsOnEven(rna, partner);
                }
                else
                {
//                  printf("C \n");
                    showGoodNucsOnOdd(rna, partner);
                }
            }
        }
    }

}


/** 
For each RNA, merge the good RNAs from all different interactions 
**/
void merge(int num)
{
#if TEST == 1
    printf("merge() \n");
#endif
    int i,j,k, currRNA, rna, partner;
    
#if TEST == 1
    for(rna=0;rna<num;rna++)
    {
        for(partner=0;partner<num;partner++)
        {
            printf("%d ", interactsWith[rna][partner]);
        }
        printf("\n");
    }
#endif
    
    for(rna=0;rna<num;rna++)
    {
//      printf("\n\nRNA = %d \n", rna);
        
        int evenOdd = data2->rnaSequences[rna]->type;
                
                
//      printf("%d interacts with: \n", rna);
        

        int container[300];
        int cc = 0;
        int interactsCount = 0;
        for(partner=0;partner<num;partner++)
        {
            if(interactsWith[rna][partner] == 1)
            {
                interactsCount++;
                //printf("%d :- ", partner);
                for(k=0;k<MAX_RNA_SIZE;k++)
                {
                    if(compressedRNAsTempHold[rna][partner][k] != -9)
                    {
                        //printf("%d ", compressedRNAsTempHold[rna][partner][k]);
                        container[cc] = compressedRNAsTempHold[rna][partner][k];
    //printf("container[%d] = %d \n", cc ,container[cc]);
                        cc++;
                    }
                }
        //      printf("\n");
            }
            
        }
        
        if(interactsCount == 0) continue;
        //now merge all the stuff in container:
        int sz = cc;
#if TEST == 1       
        printf("\verify = \n");
        for(cc=0;cc<sz;cc++)
        {
            printf("verify container[%d] = %d \n", cc, container[cc] );
        }
#endif      
        
        int cc2;
        //This loop gets rid of the duplicates (converts them to -1)
        for(cc=0;cc<sz;cc++)
        {
            for(cc2=cc+1;cc2<sz;cc2++)
            {
                if(container[cc] == container[cc2])
                    container[cc2] = -1;
            }
        }
#if TEST == 1
        printf("thing = ");
        for(cc=0;cc<sz;cc++)
        {
            printf("%d ", container[cc] );
        }
#endif

        //now sort asc if original is asc:
        
        for(cc=0;cc<sz;cc++)
        {
            for(cc2=cc+1;cc2<sz;cc2++)
            {
                if(container[cc] > container[cc2])
                {
                    
                    int temp = container[cc];
                    container[cc] = container[cc2];
                    container[cc2] = temp;
                }
            }
        }
#if TEST == 1       
        printf("\nsorted: ");
        for(cc=0;cc<sz;cc++)
            printf("%d ", container[cc] );
#endif  
        //now fill up the compressed array
        cc = 0;
        int cc3 = 0;
        //printf("\nremove -1s : \n");
        for(cc=0;cc<sz;cc++)
        {
            //printf("container[%d] = %d \t cc3 = %d\n", cc, container[cc], cc3 );
            if(container[cc] == -1)
                continue;
            data2->rnaSequences[rna]->compressedRNA[cc3+1] = container[cc];
            
            cc3++;
        }
        data2->rnaSequences[rna]->compressedRNA[0] = -5;
        
        int compLen = cc3;
#if TEST == 2   
        printf("\ncompressed array: ");
        for(cc=0;cc<=compLen;cc++)
        {
            printf("%d ", data2->rnaSequences[rna]->compressedRNA[cc]);
        }
        printf("\n");
#endif
        data2->rnaSequences[rna]->compressedLength = compLen;
        //printf("compLen = %d\n", compLen);

        //fill in expandedRNAmap
        
        int rnaLen = data2->rnaSequences[rna]->originalLength, put;
        
        for(j=0;j<=rnaLen;j++)
        {
            data2->rnaSequences[rna]->expandedRNAmap[j] = -1;
        }
//printf("sz = %d \n", sz);
//      for(j=1;j<=sz;j++)
        for(j=1;j<=compLen;j++)
        {
            //printf("j:%d, data2->rnaSequences[rna]->compressedRNA[j]-1 = %d \n", j, data2->rnaSequences[rna]->compressedRNA[j]-1);
            data2->rnaSequences[rna]->expandedRNAmap[data2->rnaSequences[rna]->compressedRNA[j]-1] = j;
        }

//if(breaker() == 1) exit(0);

        put = 0;
        for(j=0;j<=rnaLen;j++)
        {
            if(data2->rnaSequences[rna]->expandedRNAmap[j] != -1)
                put = data2->rnaSequences[rna]->expandedRNAmap[j];
            else
                data2->rnaSequences[rna]->expandedRNAmap[j] = put+1;
        }

        for(j=rnaLen;j>0;j--)
        {
            data2->rnaSequences[rna]->expandedRNAmap[j] = data2->rnaSequences[rna]->expandedRNAmap[j-1];
        }

        data2->rnaSequences[rna]->expandedRNAmap[0] = -9;
    #if TEST == 2
        printf("expandedRNAmap: \n");
        for(j=0;j<=rnaLen;j++)
        {
            printf("j: %d -> %d ", j, data2->rnaSequences[rna]->expandedRNAmap[j]);
        }
        printf("\n");
        printf("\n");
#endif
        
    }   

}







//////////////////////////////////////
////////////  New Method  ////////////
//////////////////////////////////////


//Both willl from end now

void showGoodNucsOnEven(int even, int odd)
{
#if TEST == 1
    printf("even = %d, odd = %d \n", even, odd);
    printf("C string 1 len: %d \n", data2->rnaSequences[0]->originalLength);
    printf("C string 2 len: %d \n", data2->rnaSequences[1]->originalLength);
#endif  
    
    interactsWith[even][odd] = 1;
//  printf("setting interactsWith[%d][%d] = 1 \n", even, odd);
    
//  printf("showGoodNucsOnEven: even = %d, odd = %d \n", even, odd);
    compLen = 0;
    int loc = data2->interLocs[even][odd];
    
    double **** inter = data2->rnaCollections[loc];
    
    int i,j,w,w2;
    int n1 = data2->rnaSequences[even]->originalLength;
    int n2 = data2->rnaSequences[odd]->originalLength;

    
    //printf("s1 = %s, s2 = %s \n", data2->rnaSequences[even]->string, data2->rnaSequences[odd]->string);
    //printf("s1 = %d, s2 = %d \n", data2->rnaSequences[even]->origId, data2->rnaSequences[odd]->origId);
    //printf("n1 = %d, n2 = %d \n", n1, n2);
    for(i=1;i<=n1;i++)
    {
        double min = 999;
        for(j=1;j<=n2;j++)
        {
            for(w=1;w<=25;w++)
            {
                for(w2=1;w2<=25;w2++)
                {
                    if(min > inter[i][j][w][w2])
                        min = inter[i][j][w][w2];
                }
            }
        }

        //if(min < 0)
        if(min < 9999)
        {
            ++compLen;
            //printf("%d\t%d\t%f\n", i, compLen, min);  
            compressedRNAsTempHold[even][odd][compLen] = i;
        }

    }
    
    //printf("------\n");
}


void showGoodNucsOnOdd(int odd, int even)
{
    interactsWith[odd][even] = 1;
/*  printf("setting interactsWith[%d][%d] = 1 \n", odd, even);
    
    printf("showGoodNucsOnOdd: odd = %d, even = %d \n", odd, even); */
    compLen = 0;
    int loc = data2->interLocs[odd][even];
//  printf("loc = %d\n", loc);
    double **** inter = data2->rnaCollections[loc];
//printf("here1 \n")    ;
    int i,j,w,w2;
    int n1 = data2->rnaSequences[odd]->originalLength;
//printf("here2 \n")    ;
    int n2 = data2->rnaSequences[even]->originalLength;
//printf("here3 \n")    ;
    for(i=1;i<=n1;i++)
    {
//printf("here4 \n")    ;
        double min = 999;
        for(j=1;j<=n2;j++)
        {
//printf("here 5 \n")   ;

            for(w=1;w<=25;w++)
            {
//printf("here 6  -   =%d, j=%d, w=%d  \n", i,j,w)  ;
                for(w2=1;w2<=25;w2++)
                {
                    if(min > inter[j][i][w][w2])
                        min = inter[j][i][w][w2];
                }
            }
//printf("here7 \n")    ;
        }
//printf("here8 \n")    ;
        //if(min < 0)
        if(min < 9999)
        {
            ++compLen;
//          printf("%d\t%d\t%f\n", i, compLen, min);    
            compressedRNAsTempHold[odd][even][compLen] = i;
        }
//printf("here 9 \n")   ;
    }

    //printf("------\n");
}

};

#endif