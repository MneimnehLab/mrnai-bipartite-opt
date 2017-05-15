#ifndef _RNAPROPS_H
#define _RNAPROPS_H

#include <string>
#include <vector>
#include "rnaseq.h"
#include "Weights.h"
using namespace std;


class RNAProperties
{
private:
    Weights * weights;
    int numEven;
    int numOdd;
    vector<OrigRNASeq> * origRNASequences;

public:
    RNAProperties(Weights * weights, vector<OrigRNASeq> * v);

    Weights * getWeights();
    int getNumEven();
    int getNumOdd();

    vector<int> evenSeq, oddSeq;
    vector<int> rnaLengths;

};

#endif
