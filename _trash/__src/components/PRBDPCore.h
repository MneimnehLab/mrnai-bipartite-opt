#ifndef _PRBDPCORE_H
#define _PRBDPCORE_H

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include <float.h> 
#include <cstdarg>
#include <fstream>
#include "rnaseq.h"
#include "Window.h"
#include "Config.h"
#include "../utils/MultiIndexArray.h"

class PRBDPCore
{
private:
    ofstream matrixStream;
public:

    long iterCount ;
    long H_arrayHits;
    long W_arrayHits;
    long wType;
    FILE * testfp;

    bool wrapAround;

    void backtrackNR(int * indices);
    void deallocAll();
    void reverse(char * source, char * dest);

    // extern SubEnsemble * getEnsembleFromCache(char * stringRep);
    // extern void addSubensembleToCache(SubEnsemble * subensemble);

    int gap = 0;
    int interLocs[10][10];  //interLocs[i][j] = y contains location y of interaction (i.e., rnaCollections[y])

    SingleRunConfig * data;

    vector<int> treeNodeLens;       //contains lengths of each RNA (index by RNA num)
    vector<vector<int> > elemToRNAmap;  //elemToRNAmap[i][j] contains num of the RNA at level i in element j where 1 <= j <= sum(n)
    vector<int> cutoffs;            //cutoffs[i] contains cumulative length of all RNAs prior to i.

    // int numOfRNA;        //number of RNAs
    int numOfLevels;    //number of levels
    int totalDimSize;   //product of lengths of all LEVELS
    vector<int> realLengths;    //lengths of LEVELS
    vector<int> lengths;        //lengths of LEVELS + 1
    vector<int> bases;  // bases to calculate index ( = lengths of LEVELS + 1)
    vector<int> rnaLengths; //lengths of LEVELS
    int totalLevels;

    //The following 3 arrays are linearized. indexed by l1,...,lk
    MultiIndexArray<double> H;      //The MAIN ARRAY. H.
    MultiIndexArray<short> chosenByArr; //contains "chosen by" for backtracking
    MultiIndexArray<int> chosenW;       //selected window size for backtracking
    MultiIndexArray<int> chosenW2;      //selected window size rna2 for backtracking

    // vector<double> HARRAY;
    // double * HARRAY;

    // SubEnsemble * subensemble;

    Config config;



    /// New functions

    PRBDPCore(SingleRunConfig * , int gapSize, int wType, bool wrapAround);
    ~PRBDPCore();
    void initMetaTrees(int num);
    void initLinear(int num);
    void allForLoops();
    void atStepH(int i, int * indices);



};



#endif
