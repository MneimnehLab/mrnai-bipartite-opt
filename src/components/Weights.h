#ifndef _WEIGHTS_H
#define _WEIGHTS_H

#include <string>
#include <vector>
#include "rnaseq.h"
using namespace std;


class Weights
{
private:
    // double ***** rnaupCollections;
    vector<double ****> rnaupCollections;
    vector<vector<int>> matchingMatrix;

    string dirPath;
    vector<OrigRNASeq> * origRNASequences;

    int numEven;
    int numOdd;

public:
    Weights(string dirPath, vector<OrigRNASeq> * origRNASequences);
    void Read();
    void ReadOneFile(int evenId, int oddId, int location);

    int MatchingMatrixVal(int, int);
    double **** GetWeightsTable(int);

    int getNumEven();
    int getNumOdd();

    double getWeight(int even, int odd, int evenPos, int oddPos, int evenW, int oddW);
};

#endif
