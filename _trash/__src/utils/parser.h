#ifndef _PARSE_H
#define _PARSE_H

#include <vector>
#include "../components/rnaseq.h"
using namespace std;


class Parser
{
public:
    static vector<string> ReadAndTokenizeAmpersand();
    static void ReadEvenOdd(vector<OrigRNASeq> & origRNASequences);
    static vector<OrigRNASeq> GetAndParseInput();
    static CmdLineArgs * ParseArgs(int argc, char ** argv);
};

#endif
