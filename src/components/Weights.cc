#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include "Weights.h"
#include "../utils/Arrays.h"
using std::string;

Weights::Weights(string d, vector<OrigRNASeq> * v) : dirPath(d), 
                                                    origRNASequences(v)
{

}

// 5/14: Modifying this do that we read and store all weight files
// without assuming anything about even/odd
void Weights::Read()
{
    std::cout << "Read weights from file.... \n" << std::endl;

    int numOfRNA = origRNASequences->size();
    
    // initialize matrix
    matchingMatrix.resize(numOfRNA, vector<int>(numOfRNA,-9));
    
    // do RNAup for each rna pair (we don't know about/care about even odd right now)
    int num_pairs = numOfRNA * numOfRNA;
    rnaupCollections.resize(num_pairs);
    
    int i = 0;
    for(int even=0; even<numOfRNA; even++)
    {
        for(int odd=0; odd<numOfRNA; odd++)
        {
            if(even == odd)
                continue;

            matchingMatrix[even][odd] = i;

            ReadOneFile(even, odd, i);

            i++;
        }
    }
}


void Weights::ReadOneFile(int evenId, int oddId, int location)
{
    string evenSeq = (origRNASequences->at(evenId)).string;
    string oddSeq  = (origRNASequences->at(oddId)).string;

    double ****rnaWins = Arrays::make4DTable<double>(evenSeq.length()+1, oddSeq.length()+1, 26, 26);
    
    std::stringstream fileNameStream("");
    fileNameStream << dirPath << "default-" << evenId << "_" << oddId << "_itemized.out";
    string fileName = fileNameStream.str();

    std::ifstream rnaupEnergies(fileName.c_str());

    if(rnaupEnergies.good())
    {
        string line;

        while(!rnaupEnergies.eof())
        {
            getline(rnaupEnergies, line);

            // The line may have many values; however, we only care about the first 5
            std::stringstream lineParser(line);
            
            int i1,j1,i2,j2;
            double energy;
            lineParser >> i1 >> j1 >> i2 >> j2 >> energy; 
            int w1 = j1-i1;
            int w2 = j2-i2;

            rnaWins[j1][j2][w1][w2] = energy;
        }

        std::cout << "Read file " << fileName.c_str() << endl;
        
        rnaupCollections[location] = rnaWins;
    }
    else
    {
        std::cerr << "Cannot read input file " << fileName << endl;
        exit(0);
    }
}


int Weights::MatchingMatrixVal(int i, int j)
{
    return matchingMatrix[i][j];
}

double **** Weights::GetWeightsTable(int i)
{
    return rnaupCollections[i];
}

int Weights::getNumEven()
{
    return numEven;
}

int Weights::getNumOdd()
{
    return numOdd;
}

double Weights::getWeight(int even, int odd, int evenPos, int oddPos, int evenW, int oddW)
{
    int loc = matchingMatrix[even][odd];
    return rnaupCollections[loc][evenPos][oddPos][evenW][oddW];
}
