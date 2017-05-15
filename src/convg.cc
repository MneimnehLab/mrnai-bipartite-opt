#include <algorithm>
#include <iostream>
#include <list>
#include <numeric>
#include <random>
#include <vector>
#include <exception>

#include "components/rnaseq.h"
#include "components/Weights.h"
#include "components/Config.h"
#include "components/PRBDPCore.h"
#include "utils/parser.h"


void relabelSeqTypes(vector<OrigRNASeq> & rnaSequences, vector<int> ordering);
vector<vector<int> > getNeighbors(const vector<int> & ordering);
vector<int> randomEvenOddVector(int num);


int main(int argc, char *argv[])
{
    
    // CmdLineArgs * args = Parser::ParseArgs(argc, argv);
    vector<OrigRNASeq> origRNASequences = Parser::GetAndParseInput();
    int numOfRNA = origRNASequences.size();

    string dir = "output/";
    Weights weights(dir, &origRNASequences);
    weights.Read();

    // vector<int> ordering {0,-1,1,2,3};
    vector<int> ordering = randomEvenOddVector(numOfRNA);
    vector<int> minOrder;
    double minEnergy;
    Config minEnergyConfig;

    bool done = false;

    while(!done)
    {
        bool foundBetterStructure = false;
        for(auto newOrdering : getNeighbors(ordering))
        {
            cout << "New ordering: ";
            for(auto e : newOrdering)
                cout << e << " ";
            cout << endl;

            relabelSeqTypes(origRNASequences, newOrdering);
            
            RNAProperties props(&weights, &origRNASequences);

            PRBDPCore dpAlgo(&props); 

            double energy = dpAlgo.getMinEnergy();

            if(energy < minEnergy)
            {
                minEnergy = energy;
                minOrder = newOrdering;
                foundBetterStructure = true;
                minEnergyConfig = dpAlgo.getResultConfig();
            }
            cout << endl;
        }

        if(!foundBetterStructure)
            done = true;

        cout << endl << endl;
    }

    cout << endl << endl;

    cout << "Min Energy = " << minEnergy << endl;
    cout << "Min Energy Structure = " << endl;
    cout << minEnergyConfig << endl;
    cout << "Min ordering: ";
            for(auto e : minOrder)
                cout << e << " ";
            cout << endl;

    cout << endl;

    return 0;
}



// Expect a -1 in the ordering vector, which is the barrier
void relabelSeqTypes(vector<OrigRNASeq> & rnaSequences, vector<int> ordering)
{
    vector<OrigRNASeq> temp;
    int type = 0;
    for(int id : ordering)
    {
        if(id == -1)    // barrier
        {
            // flip type
            type = 1;
            continue;
        }
        rnaSequences[id].type = type;
        // temp.push_back(rnaSequences[id]);
    }

    if(type == 0)
        throw std::logic_error("No barrier in ordering!\n");
}

vector<vector<int> > getNeighbors(const vector<int> & ordering)
{
    vector<vector<int> > neighbors;
    int N = ordering.size();
    int bar_idx = -1;
    for(int i=0; i<N; i++)
        if(ordering[i] == -1)
            bar_idx = i;

    if(bar_idx == -1)
        throw std::logic_error("No barrier in ordering!\n");

    auto itI=ordering.begin();
    for(int i=0 ; i<bar_idx + 1; itI++, i++)    // remove + 1 if we don't want to include the barrier
    {
        auto itJ=ordering.begin() + bar_idx ;  // add + 1 if we don't want to include the barrier
        for(int j=bar_idx; j<N; itJ++, j++)     // add + 1 to j if we don't want to include the barrier
        {
            if(i == j) continue;

            vector<int> nbr;
            // insert mid section and reverse
            nbr.insert(nbr.end(), itI, itJ+1);
            std::reverse(nbr.begin(), nbr.end());
            
            // append front and back
            nbr.insert(nbr.begin(), ordering.begin(), itI);
            nbr.insert(nbr.end(), itJ+1, ordering.end());

            if(*nbr.begin() == -1 or *(nbr.end()-1) == -1)
                continue;

            neighbors.push_back(nbr);
        }
    }
    return neighbors;
}

vector<int> randomEvenOddVector(int num)
{
    vector<int> ordering(num);
    std::iota(ordering.begin(), ordering.end(), 0);
    std::shuffle(ordering.begin(), ordering.end(), std::mt19937{std::random_device{}()});

    // insert barrier in some non-endpoint position
    auto x = std::mt19937{std::random_device{}()};
    int pos = x() % (num-1) + 1;
    ordering.insert(ordering.begin() + pos, -1);

    return ordering;
}

