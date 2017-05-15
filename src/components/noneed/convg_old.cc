#include <iostream>
#include "components/rnaseq.h"
#include "components/Weights.h"
#include "components/Cache.h"
#include "components/Chain.h"
#include "components/Config.h"
#include "components/PRBDPCore.h"
#include "components/LevelGroupProcessor.h"
#include "utils/parser.h"

#define PROG_TYPE "findone_"

void printMilliSecs(){}

void printTime() {}

int main(int argc, char *argv[])
{
    
    CmdLineArgs * args = Parser::ParseArgs(argc, argv);
    vector<OrigRNASeq> origRNASequences = Parser::GetAndParseInput();
    int numOfRNA = origRNASequences.size();

    // string dir = "../rnaup_weights/output/";
    // string dir = "../output/";
    string dir = "output/";
    Weights weights(dir, &origRNASequences);
    weights.Read();

    RNAProperties props(&weights, &origRNASequences);

    PRBDPCore dpAlgo(&props); 
    
















    /*
    Chain chain(&origRNASequences, numOfRNA);
    // vector<int>v {0,1,3,2};
    // Chain chain = Chain::makeGivenChain(&origRNASequences, v);
    chain.printFlatStruct();
    chain.determineStruct();

    if(args->k > chain.getBBHeight())
        args->k = chain.getBBHeight();
    LevelGroupProcessor lgProc(args->k, &chain, &weights, args->loopAround);
    
    vector<vector<int> > x = lgProc.GetAllGroupings();

    double minEnergy = 0;
    Config minConfig;

    for(int groupNum=0; groupNum<x.size(); groupNum++)
    {
        cout << "Computing Opt for Grouping " << (groupNum+1) 
             << "/" << x.size() << endl
             << "------------------------------------" << endl;
        // cout << "\n\n\nin group " << groupNum << endl << endl;
        Config c;
        lgProc.prepareWithSubsets(groupNum, c);
        if(minEnergy > c.totalWeight())
        {
            minEnergy = c.totalWeight();
            minConfig = c;
        }

        // cout << "\n\n\ndone with group " << groupNum << endl << endl;

        cout << endl << endl;
    }

    cout << endl
         << "Minimum Energy: " << minEnergy << endl
         << "Structure: " << endl;
    for(auto w: minConfig)
        w.prettyPrint();
    */

    return 0;
}



