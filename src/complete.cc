#include <iostream>
#include "components/rnaseq.h"
#include "components/Weights.h"
#include "components/Config.h"
#include "components/PRBDPCore.h"
#include "utils/parser.h"




int main(int argc, char *argv[])
{
    
    // CmdLineArgs * args = Parser::ParseArgs(argc, argv);
    vector<OrigRNASeq> origRNASequences = Parser::GetAndParseInput();
    Parser::ReadEvenOdd(origRNASequences);
    
    
    string dir = "output/";
    Weights weights(dir, &origRNASequences);
    weights.Read();

    RNAProperties props(&weights, &origRNASequences);

    PRBDPCore dpAlgo(&props); 

    double energy = dpAlgo.getMinEnergy();
    Config config = dpAlgo.getResultConfig();

    cout << "Min Energy = " << energy << endl;
    cout << "Min Energy Structure = " << endl;
    cout << config << endl;
    
    cout << endl;

    return 0;
}