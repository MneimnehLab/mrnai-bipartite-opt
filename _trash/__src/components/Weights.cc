#include <iostream>
#include <sstream>
#include <fstream>
#include "Weights.h"
#include "../utils/Arrays.h"
using std::string;

Weights::Weights(string d, vector<OrigRNASeq> * v) : dirPath(d), 
                                                    origRNASequences(v)
{

}

void Weights::Read()
{
    std::cout << "Read weights from file.... \n" << std::endl;

    int numOfRNA = origRNASequences->size();
    // int xi;

    int i;
    short evenSeq[10], e = 0;
    short oddSeq[10], o = 0;

    // initialize matrix
    for(i=0;i<100;i++)
        matchingMatrix[i/10][i%10] = -9;    

    // find even or odd
    for(i=0;i<numOfRNA;i++)
    {
        if((origRNASequences->at(i)).type == 0)
            evenSeq[e++] = i;
        else if((origRNASequences->at(i)).type == 1)
            oddSeq[o++] = i;    
    }

    //do RNAup for each even/odd pair
    int pairs = o * e;

    rnaupCollections = (double*****)malloc(sizeof(double ****) * pairs);

    i = 0;
    for(int oi=0;oi<o;oi++)
    {
        for(int ei=0;ei<e;ei++)
        {

            int even = evenSeq[ei];
            int odd  = oddSeq[oi];

            matchingMatrix[even][odd] = i;
            matchingMatrix[odd][even] = i;

            // int evenLen = (origRNASequences->at(even)).string.length();
            // int oddLen = (origRNASequences->at(odd)).string.length();
            
            ReadOneFile((origRNASequences->at(even)).origId, (origRNASequences->at(odd)).origId, i);
            
            // if(wType == 1)
            //  doSubAddFor2RNAs_sameW(i, evenLen, oddLen, even, odd);
            // else
            //  doSubAddFor2RNAs_twoW(i, evenLen, oddLen, even, odd);

            i++;
        }
    }
}

// void readPirna(int evenId, int oddId, int n1, int n2, int location, double ***** rnaupCollections, char * fileName2)
void Weights::ReadOneFile(int evenId, int oddId, int location)
{
    string evenSeq = (origRNASequences->at(evenId)).string;
    string oddSeq  = (origRNASequences->at(oddId)).string;

    double ****rnaWins = Arrays::make4DTable<double>(evenSeq.length()+1, oddSeq.length()+1, 26, 26);
    
    std::stringstream fileNameStream("");
    fileNameStream << dirPath << "default-" << evenId << "_" << oddId << "_itemized.out";

    string fileName = fileNameStream.str();
    // sprintf(fileName, "../rnaup_weights/output/default-%d_%d.weights", evenId, oddId);
    
    // int count = 0;
    // std::FILE * file = fopen(fileName.c_str(), "rt");

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

        printf("Read file %s\n", fileName.c_str());
        rnaupCollections[location] = rnaWins;
        
    }
    /*
    if ( file != NULL )
    {
        char line [ 500 ]; 
        while ( fgets ( line, sizeof line, file ) != NULL ) 
        {       
            
            int i1,j1,i2,j2;
            double energy;
            //printf("%s\n", line);
            sscanf (line, "%d\t%d\t%d\t%d\t%lf\n", &i1,&j1,&i2,&j2,&energy);
            
            int w1 = j1-i1;
            int w2 = j2-i2;
            rnaWins[j1][j2][w1][w2] = energy;
            
            // std::printf("%d, %d, %d, %d = %f \n", j1,j2,w1,w2, energy);
            //printf("%d, %d, %d, %d = %f \n", i1,j1,i2,j2, energy);
            
        }
        fclose ( file );
        printf("Read file %s\n", fileName.c_str());
        rnaupCollections[location] = rnaWins;
    }
    */
    else
    {
        printf("File read error! Could not read file %s \n", fileName.c_str());
        exit(1);
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
