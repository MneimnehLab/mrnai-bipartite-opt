#ifndef _WEIGHTS_H
#define _WEIGHTS_H

#include <string>
#include <vector>
#include "rnaseq.h"
using namespace std;


class Weights
{
private:
    double ***** rnaupCollections;
    int matchingMatrix[10][10];

    string dirPath;
    vector<OrigRNASeq> * origRNASequences;

public:
    Weights(string dirPath, vector<OrigRNASeq> * origRNASequences);
    void Read();
    void ReadOneFile(int evenId, int oddId, int location);

    int MatchingMatrixVal(int, int);
    double **** GetWeightsTable(int);
};

#endif
