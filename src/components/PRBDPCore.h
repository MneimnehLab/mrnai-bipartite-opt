#ifndef _PRBDPCORE_H
#define _PRBDPCORE_H

#include "RNAProperties.h"
#include "Window.h"
#include "Weights.h"
#include "Config.h"
#include "../utils/MultiIndexArray.h"

class PRBDPCore
{
private:
    struct RNA_Pair
    {
        int even, odd;
    };

    Weights * weights;
    int numEven;
    int numOdd;

    vector<int> evenRNAs, oddRNAs;

    Config result_config;

    double minEnergy;
    
    // int gap = 0;

    int totalDimSize;       //product of lengths of all LEVELS
    vector<int> lengths;    //lengths of RNAs + 1
    int totalLevels;


    //The following 3 arrays are linearized. indexed by l1,...,lk
    MultiIndexArray<double> H;              //The MAIN ARRAY. H.
    MultiIndexArray<RNA_Pair> chosenByArr;  //contains "chosen by" for backtracking
    MultiIndexArray<int> chosenW;           //selected window size for backtracking
    MultiIndexArray<int> chosenW2;          //selected window size rna2 for backtracking

    void backtrackNR(int * indices);
    void initLinear(int num);
    void allForLoops();
    void atStepH(int i, int * indices);

public:
    
    Config getResultConfig();

    double getMinEnergy();
    
    PRBDPCore(RNAProperties *rnaProps);
    ~PRBDPCore();
    
};



#endif
