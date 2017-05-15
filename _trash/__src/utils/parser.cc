#include <sstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include "parser.h"
#include "../components/rnaseq.h"

using std::vector;
using std::string;
using std::cout;
using std::endl;


// following two functions thanks to "iain"! http://stackoverflow.com/a/868894/554658

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}


vector<string> Parser::ReadAndTokenizeAmpersand()
{
    // This function can be used for RNA sequences as well as names
    // both types are separated by '&'s
    
    string input;
    std::getline(std::cin, input);
    
    vector<string> sequences;
    std::stringstream parserStream(input);
    for(string sequence; 
        getline(parserStream, sequence, '&'); 
        sequences.push_back(sequence) );

    return sequences;
}

void Parser::ReadEvenOdd(vector<OrigRNASeq> & origRNASequences)
{
    std::cout << "Enter even1,even2,..|odd1,odd2,...: ";

    string line;
    std::getline(std::cin, line);
    
    int readingType = 0;
    
    for(int i=0; i<line.length(); i++)
    {
        switch(line[i])
        {
            case '|':
                readingType = 1;
                break;
            case ',':
                break;
            default:
                origRNASequences[line[i]-48].type = readingType;
        }
    }
}

vector<OrigRNASeq> Parser::GetAndParseInput()
{
    
    std::cout << "Enter Sequences: ";
    vector<string> sequences = Parser::ReadAndTokenizeAmpersand();

    std::cout << "Enter Names: ";
    vector<string> names = Parser::ReadAndTokenizeAmpersand();

    vector<OrigRNASeq> origRNASequences;
    for(int i=0; i<sequences.size(); i++)
    {
        string seq = sequences[i];
        string name = names[i];

        OrigRNASeq o = {i, seq, (int)seq.length(), name};
        origRNASequences.push_back( o );
    }

    Parser::ReadEvenOdd(origRNASequences);

    return origRNASequences;

}


CmdLineArgs * Parser::ParseArgs(int argc, char ** argv)
{
    if(cmdOptionExists(argv, argv+argc, "-h"))
    {
        cout << "Usage: ./findone  [-p ] flag: parallelizes processing" << endl;
        cout << "                  [-k   number of levels (1/epsilon)]" << endl;
        cout << "                  [-t   number of trials]" << endl;
        cout << "                  [-w   windows sizes(1=same or 2=diff)]" << endl;
        cout << "                  [-gb  fill level gaps (0=no or 1=yes)]" << endl;
        cout << "                  [-gs  window gap size]" << endl;
        cout << "                  [-gu  allow GU pairs (0=no or 1=yes)]" << endl;
        cout << "                  [-f   input file name]" << endl;
        cout << "                  [-r   rnaup out file name]" << endl;
        cout << "                  [-nl   (no loop around)]" << endl;

        exit(1);
    }

    CmdLineArgs * args = new CmdLineArgs;
    
    int i;
    int k, t, gbool, gsize, w, gu;
    for (i=1; i<argc; i++) 
    {
        if (argv[i][0]=='-') 
            switch ( argv[i][1] )
            {
                case 'h':
                    cout << "Usage: ./findone  [-p ] flag: parallelizes processing" << endl;
                    cout << "                  [-k   number of levels (1/epsilon)]" << endl;
                    cout << "                  [-t   number of trials]" << endl;
                    cout << "                  [-w   windows sizes(1=same or 2=diff)]" << endl;
                    cout << "                  [-gb  fill level gaps (0=no or 1=yes)]" << endl;
                    cout << "                  [-gs  window gap size]" << endl;
                    cout << "                  [-gu  allow GU pairs (0=no or 1=yes)]" << endl;
                    cout << "                  [-f   input file name]" << endl;
                    cout << "                  [-r   rnaup out file name]" << endl;
                    cout << "                  [-nl ] (no loop around)]" << endl;

                    exit(1);
                //parallel
                case 'p':  
                    args->parallel = 1;
                    break;

                //k
                case 'k':  
                    k = atoi(argv[i+1]);
                    if(k == 0) {cerr << "k has to be non zero \n" << endl; exit(-1); }
                    args->k = k ;
                    break;
    
                //number of trials
                case 't':  
                    t = atoi(argv[i+1]);
                    if(t == 0) {cerr << "t has to be non zero \n" << endl; exit(-1); }
                    args->trials = t ;
                    break;
                
                case 'w':  
                    w = atoi(argv[i+1]);
                    if(w < 1 || w > 2) {cerr << "w has to be 1 or 2 \n" << endl; exit(-1); }
                    args->winSize = w ;
                    break;
                
                case 'g':  
                    
                    if(argv[i][2] == 'b')
                    {
                        gbool = atoi(argv[i+1]);
                        args->fillGaps = gbool;
                    }
                    else if(argv[i][2] == 's')
                    {
                        gsize = atoi(argv[i+1]);
                        args->gapSize = gsize;
                    }
                    else if(argv[i][2] == 'u')
                    {
                        gu = atoi(argv[i+1]);
                        args->GU = gu;
                    }
                    break;

                case 'f':
                    // printf("here: %s \n", argv[i+1]);
                    args->fileName = string(argv[i+1]);
                
                case 'r':
                    printf("rnaup out: %s \n", argv[i+1]);
                    args->rnaupOut = string(argv[i+1]);
                
                default: 
                    ;
            } 
    }

    if(cmdOptionExists(argv, argv+argc, "-nl"))
    {
        args->loopAround = false;   
    }

    cerr << "Using k " << args->k;
    cerr << "Using parallel =  " << args->parallel << endl;
    cerr << "using trials =  " << args->trials << endl;
    cerr << "using gaps =  " << args->fillGaps << endl;

    return args;
}





/*
CmdLineArgs * Parser::parseArgs(int argc, char ** argv)
{
    if(cmdOptionExists(argv, argv+argc, "-h"))
    {
        cout << "Usage: ./main  [-f output_filename]";
        cout << "               [-s number of samples]";

        exit(1);
    }

    CmdLineArgs * args = new CmdLineArgs;   // Constructor will set default args

    if(cmdOptionExists(argv, argv+argc, "-p"))
        args->parallel = 1;
    
    char * cmdChr = getCmdOption(argv, argv + argc, "-k");
    if (cmdChr)
        args->k = stoi(string(cmdChr));

    cmdChr = getCmdOption(argv, argv + argc, "-k");
    if (cmdChr)
        args->k = stoi(string(cmdChr));
    if(args->k <= 0) {cerr << "k has to be positive \n" << endl; exit(-1); }

    cmdChr = getCmdOption(argv, argv + argc, "-t");
    if (cmdChr)
        args->trials = stoi(string(cmdChr));
    if(args->trials <= 0) {cerr << "k has to be positive \n" << endl; exit(-1); }

    cmdChr = getCmdOption(argv, argv + argc, "-w");
    if (cmdChr)
        args->w = stoi(string(cmdChr));
    if(args->w != and args->w != 2) {cerr << "w has to be 1 or 2 \n" << endl; exit(-1); }

    cmdChr = getCmdOption(argv, argv + argc, "-gb");
    if (cmdChr)
        args->fillGaps = stoi(string(cmdChr));
    if(args->fillGaps != 0 and args->fillGaps != 1) {cerr << "-gb has to be 0 or 1\n" << endl; exit(-1); }


    
}
*/


// int main()
// {
//  Parser::GetAndParseInput();

//  return 0;
// }
