#include <iostream>
#include "components/rnaseq.h"
#include "components/Weights.h"
#include "components/Cache.h"
#include "components/Chain.h"
#include "components/Config.h"
#include "components/LevelGroupProcessor.h"
#include "utils/parser.h"

#define PROG_TYPE "findone_"

void printMilliSecs(){}

void printTime() {}

void perms(vector<vector<int>> &allPerms, vector<int> elems, vector<int> collected = vector<int>() )
{
    if(elems.size() == 0)
    {
        allPerms.push_back(collected);
        return;
    }

    for(auto newFirst : elems)
    {
        vector<int> remaining;
        
        for(auto ej : elems)
            if(ej != newFirst)
                remaining.push_back(ej);
        
        vector<int> collected2(collected);
        collected2.push_back(newFirst);
        
        perms(allPerms, remaining, collected2);
    }
}


int main(int argc, char *argv[])
{
    
    CmdLineArgs * args = Parser::ParseArgs(argc, argv);
    vector<OrigRNASeq> origRNASequences = Parser::GetAndParseInput();
    int numOfRNA = origRNASequences.size();

    Cache matchingCache;

    // string dir = "../rnaup_weights/output/";
    // string dir = "../output/";
    string dir = "output/";
    Weights weights(dir, &origRNASequences);
    weights.Read();

    vector<int> v;
    for(int i=0; i< origRNASequences.size(); i++)
        v.push_back(i);

    vector<vector<int>> allPerms;
    perms(allPerms, v);
    
    double minEnergy = 0;
    Config minConfig;

    for(auto permutation : allPerms)
    {
        int k = args->k;

        cout << "Permutation: ";
        for(auto e: permutation)
            cout << e << " " ;
        cout << endl;
        
        Chain chain = Chain::makeGivenChain(&origRNASequences, permutation);
        chain.printFlatStruct();
        chain.determineStruct();

        if(k > chain.getBBHeight())
            k = chain.getBBHeight();
        
        LevelGroupProcessor lgProc(k, &chain, &weights, args->loopAround);
        
        vector<vector<int> > x = lgProc.GetAllGroupings();
        for(int groupNum=0; groupNum<x.size(); groupNum++)
        {
            Config c;
            lgProc.prepareWithSubsets(groupNum, c);
            if(minEnergy > c.totalWeight())
            {
                minEnergy = c.totalWeight();
                minConfig = c;
            }
        }

        cout << endl << endl;
    }


    cout << endl
         << "Minimum Energy: " << minEnergy << endl
         << "Structure: " << endl;
    for(auto w: minConfig)
        w.prettyPrint();
    
    return 0;
}



