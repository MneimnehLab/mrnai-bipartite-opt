#ifndef _RNASEQ_H
#define _RNASEQ_H

#include <string>
#include <vector>
using namespace std;

#define MAX_RNA_SIZE 250
#define MAX_NUM_RNA 10
#define GAP 4

typedef struct OrigRNASeq
{
    int origId;
    string string;
    int originalLength;
    std::string name;
    int type;
    
    
} OrigRNASeq;


typedef struct CmdLineArgs
{
    int  parallel;
    int  k;
    int  fillGaps;
    int  trials;
    int  gapSize;
    int  winSize;
    int  GU;
    string rnaupOut;
    string fileName;
    bool loopAround;

    CmdLineArgs() 
    {
        // default args values:
        parallel = 0;
        k = 2;
        trials = 1;
        fillGaps = 1;
        gapSize = 4;
        winSize = 1;
        GU = 0;
        fileName = "default";
        rnaupOut = "default";
        loopAround = true;
    };
} CmdLineArgs;



typedef struct RNASeq
{
    int id;
    std::string string;
    std::string name;
    
    int * compressedRNA;
    int * expandedRNAmap;
    vector<int> interactsWith;
    
    //int * parentOf[10];
    
    int compressedLength;
    int originalLength;
    
    int type;   //even or odd
    
    int origId;
    
} RNASeq;



typedef struct InteractWeight
{
    double **** weights;
} InteractWeight;



typedef struct SingleRunConfig
{
    int numOfRNA;
    vector<string> strArray1;
    int setup[MAX_NUM_RNA][MAX_NUM_RNA];
    int parentOf[MAX_NUM_RNA][MAX_NUM_RNA];
    int interLocs[MAX_NUM_RNA][MAX_NUM_RNA];
    int totalLevels;
    // vector<double ****> rnaCollections;  // RNAs in this run
    double **** rnaCollections[50]; // RNAs in this run
    int rnaCollSize;
    // vector<int> newId2OldId;
    int newId2OldId[50];
    std::string stringRep;
    int gapSize;

    vector<RNASeq *> rnaSequences;

} SingleRunConfig;

typedef struct InterDims
{
    int n1;
    int n2;
} InterDims;






typedef struct Matching
{
    double energy;
    
    
} Matching;



/*
typedef struct OneRNAup
{
    double
} OneRNAup;
*/

typedef struct OneInteractionStructure
{
    int topSeq, botSeq, startOnTopSeq, indexOnTopSeq, startOnBotSeq, indexOnBotSeq;
    int trueTopSeq, trueBotSeq;
    double weight;
} OneInteractionStructure;

/*
typedef struct SubEnsemble
{
    // OneInteractionStructure * interactions[100];
    vector<OneInteractionStructure> interactions;
    int count;
    double subenergy;
    // char * stringRep;
    string stringRep;
    int firstRowRNAs[MAX_NUM_RNA], lastRowRNAs[MAX_NUM_RNA];
    
} SubEnsemble;

typedef struct Ensemble
{
    int k;
    int first;
    int count;
    double energy;
    SubEnsemble * subs[10];
} Ensemble;
*/


typedef struct RNAupOutput
{
    
} RNAupOutput;


#endif