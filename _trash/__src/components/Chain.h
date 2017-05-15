//Even_Odd_Barrier :- EOB
//typedef struct EOB EOB;
//typedef struct Pointer Pointer;


#ifndef _CHAIN_H
#define _CHAIN_H

#include <string>
#include <vector>
#include <vector>
#include "rnaseq.h"
using namespace std;


class Chain
{
public:
    typedef struct EOB
    {
        int type;   //0: even, 1: odd
        int origId;
        
        struct EOB * next;
    } EOB;

    // Not sure why I had this, but need to keep this to prevent
    // code from breaking
    typedef struct Pointer
    {
        EOB * elem;
    } Pointer;

public:
    Pointer * head;

    void destroyStruct(Pointer * head);
    void destroyStruct(EOB * n, int i);

    vector<OrigRNASeq>* origRNASequences;
    
private:
    Chain(Pointer *, vector<OrigRNASeq> * origRNASequences);

public:
    // Chain();
    
    Chain(vector<OrigRNASeq> * origRNASequences, int numOfRNA);
    ~Chain();

    static Pointer * createEmptyStructure();
    static Pointer * makeSomeBinStruct(vector<OrigRNASeq> * origRNASequences, int numOfRNA);
    static Chain makeGivenChain(vector<OrigRNASeq> * origRNASequences, vector<int> sequence);
    static Chain makeRandomChain(vector<OrigRNASeq> * origRNASequences);
    // 
    // static Chain makeRandomBinStruct(OrigRNASeq ** origRNASequences, int numOfRNA);

    // Modifies the chain, doesn't make a copy
    void flipSubChain(int start, int end);
    void printChain();
    void determineStruct();
    int getBBHeight();
    SingleRunConfig * generateSubsetConfigFromBins(int from, int to);
    void printFlatStruct2File(std::FILE * file);
    void printFlatStruct();


    

};

#endif
