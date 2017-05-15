#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <string>
#include <limits>
#include "PRBDPCore.h"
#include "RNAProperties.h"
#include "Window.h"
#include "Weights.h"
#include "Config.h"
#include "../utils/MultiIndexArray.h"

using namespace std;

#define TEST 0

PRBDPCore::PRBDPCore(RNAProperties * rnaProps) :
    weights(rnaProps->getWeights()),
    numEven(rnaProps->getNumEven()),
    numOdd(rnaProps->getNumOdd()),
    evenRNAs(rnaProps->evenSeq),
    oddRNAs(rnaProps->oddSeq)
{
    totalLevels = numEven + numOdd;
    lengths.resize(totalLevels);
    
    for(int i=0;i<totalLevels;i++)
        lengths[i] = rnaProps->rnaLengths[i] + 1;
        
    initLinear(totalLevels);
    allForLoops();
}
   


// All stuff related to linearization initialization
void PRBDPCore::initLinear(int num)
{
    // find out the total number of elements the grand array is going to hold
    // this is = len(level1) * ... * len(level_k)
    totalDimSize = 1;
    for(int i=0;i<totalLevels;i++)
        totalDimSize *= lengths[i];
    
    H = MultiIndexArray<double>::CreateMultiIndexArray(lengths);
    chosenByArr = MultiIndexArray<RNA_Pair>::CreateMultiIndexArray(lengths);
    chosenW = MultiIndexArray<int>::CreateMultiIndexArray(lengths);
    chosenW2 = MultiIndexArray<int>::CreateMultiIndexArray(lengths);
}



void PRBDPCore::allForLoops()
{

    auto begin = chrono::high_resolution_clock::now();  

    // recursion base cases
    H[0] = 0;   
    chosenByArr[0] = {-1,-1};

    // start i from where all indices >= 1
    int block = totalDimSize/10, percent = 0;
    int * indices = new int[totalLevels];
    
    cout << "totalDimSize = " << totalDimSize << endl;
    cout << "Completed: " << endl;
    for(int i=1;i<totalDimSize;i++)
    {
        
        H.revHash(i, indices);
        
        atStepH(i, indices);
        
        if(i % block == 0) 
        { 
            percent += 10;  
            printf("%d %% .. ", percent); 
            fflush(stdout);
        };

    }
    cout << endl;

#if TEST == 1
    printf("iterCount = %ld\n", iterCount++);
    printf("H_arrayHits = %d\n", H_arrayHits++);
    printf("W_arrayHits = %d\n", W_arrayHits++);
    
#endif
    int * finalIndices = new int[totalLevels];
    
    for(int i=0;i<totalLevels;i++)
        finalIndices[i] = lengths[i] - 1;

    minEnergy = H[totalDimSize-1];
    cout << "Energy: " << minEnergy << endl;
    
    backtrackNR(finalIndices);      

    delete[] finalIndices;
    delete[] indices;

    auto end = chrono::high_resolution_clock::now();    
    auto dur = end - begin;
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur).count();
    
    cout << "Time after all H:" << (double)ms << "ms" << endl;

}





void PRBDPCore::atStepH(int theIndex, int * indices)
{
    // Note that "theIndex" is the hash of "indices", but since
    // we will require both in this function, we pass both as args
    // instead of recomputing the hash each time

    #if TEST == 1   
        int ci;
        for(ci=0;ci<totalLevels;ci++)
            // This check is only for testing...
            // if "theIndex" is unhashed correctly, no index will be < 0
            if(indices[ci] < 0)
            {
                cout << "val less than zero" << endl;
                exit(0);
                return;
            }
        
      
        for(ci=0;ci<totalLevels;ci++)
            cout << indices[ci] << " ";
        cout << endl;
    #endif
    
    double val, min = std::numeric_limits<double>::max();
    int argMinI, argMinW, argMinW2;
    RNA_Pair chosenBy = {-1, -1};



    // First part: move one peg on each level: 
    // H(i-1,j,k,..), H(i,j-1,k,...), H(i,j,k-1,...), etc
    for(int i=0; i<totalLevels; i++)
    {
        // Change one index instead of comptuing from scratch
        indices[i]--;
        if(indices[i] < 0)
        {
            indices[i]++;
            continue;
        }
        val = H[indices];
        

        if(val < min)
        {
            min = val;
            argMinI = i;
            chosenBy = {i, -1};

        #if TEST == 1
            cout << "Found min in part 1 \n";
            for(ci=0;ci<totalLevels;ci++)
                cout << indices[ci] << " ";
            cout << "\nWith min = " << min;
            cout << "\nchosen by = " <<  chosenBy.even << endl;
        #endif
        }
        
        // restore original value so we can reuse this.
        indices[i]++;
    }

    // Second part: 
    // here we look at all possible (even, odd) pairs
    
    for(int even : evenRNAs)
        for(int odd : oddRNAs)
        {
            int eIndex = even;
            int oIndex = odd;

            int currEvenPos = indices[eIndex];
            int currOddPos =  indices[oIndex];

            for(int w1=1; w1<25; w1++)
            {
                for(int w2=1; w2<25; w2++)
                {
                    int substrEvenPos = currEvenPos - w1 - 1;
                    int substrOddPos  = currOddPos  - w2 - 1;

                    if(substrEvenPos >= 0 && substrOddPos >= 0)
                    {
                        double windowWeight = weights->
                                getWeight(even, odd,  currEvenPos, currOddPos, w1, w2);

                        indices[eIndex] = substrEvenPos;
                        indices[oIndex] = substrOddPos;

                        val = windowWeight + H[indices];

                        if(val < min)
                        {
                            min = val;
                            chosenBy = {eIndex, oIndex};
                            argMinW = w1;
                            argMinW2 = w2;
                            #if TEST == 1 
                                cout << "Found min in part 2 \n";
                                cout << "val < min = " << min << endl;
                                cout << "update chosen by to " << chosenBy.even << ", " << chosenBy.odd << endl;
                            #endif
                        }
                    }
                }
            }

        // reset indices
        indices[eIndex] = currEvenPos;
        indices[oIndex] = currOddPos;

    }

    #if TEST == 1 
        cout << "Setting H[" << theIndex << "] = " << min << endl;
        cout << "Setting chosenByArr[" << theIndex << "] = " << chosenBy.even << ", " << chosenBy.odd << endl;
    #endif

    H[theIndex] = min;
    chosenByArr[theIndex] = chosenBy;
    chosenW[theIndex] = argMinW;
    chosenW2[theIndex] = argMinW2;
}




void PRBDPCore::backtrackNR(int * indices)
{   

    for(int ci=0; ci<totalLevels; ci++)
        if(indices[ci] < 0)
            return;
    
    RNA_Pair cb = chosenByArr[indices];
    
    if(cb.even == -1)
        return; // base case initialized to {-1,-1}
    else if(cb.odd == -1)
    {
        indices[cb.even]--;
        backtrackNR(indices);
    }
    // else if we are dealing with an interaction (if cb == 0 return)
    else
    {
        int eIndex = cb.even;
        int oIndex = cb.odd;

        short w1 = chosenW[indices];
        short w2 = chosenW2[indices];

        int currEvenPos = indices[eIndex];
        int currOddPos =  indices[oIndex];

        double weight = weights-> getWeight(eIndex, oIndex,  currEvenPos, currOddPos, w1, w2);
        int start_even = currEvenPos - w1;
        int start_odd  = currOddPos  - w2;
            
        printf("Interaction b/w RNA %d & %d : [%d,%d] & [%d,%d] \n", eIndex, oIndex, start_even, currEvenPos, start_odd, currOddPos);
        
        Window win(eIndex, oIndex, currEvenPos, currOddPos, w1, w2, weight);
        
        result_config.add(win);

        indices[eIndex] -= w1 + 1;
        indices[oIndex] -= w2 + 1;
        backtrackNR(indices);
    }
}

Config PRBDPCore::getResultConfig()
{
    return result_config;
}


PRBDPCore::~PRBDPCore()
{
}

double PRBDPCore::getMinEnergy()
{
    return minEnergy;
}
