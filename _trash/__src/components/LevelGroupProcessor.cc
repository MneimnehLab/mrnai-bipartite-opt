// Here "Group" = "Window of levels"

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include "LevelGroupProcessor.h"
#include "rnaseq.h"
#include "Chain.h"
#include "PRBDPCore.h"
#include "TempCompression.h"
#include "Window.h"
#include "Config.h"

#include "../utils/Arrays.h"
using namespace std;


LevelGroupProcessor::LevelGroupProcessor(int k_, Chain* chain_, Weights* weights_, bool loopAround) :
        k(k_), chain(chain_), weights(weights_), loopAround(loopAround) { }

vector<vector<int> > LevelGroupProcessor::GetAllGroupings()
{
    int levels = chain->getBBHeight();
    
    // (member) allGroups[i] = vector containing. num of levels in each 
    //                         group when first group has k-i levels
    
    // In round 0, first group has k levels, and in each subsequent round, 
    // slide window up, so first group has k - 1 levels, then k - 2, and so on.
    // Variable "first" = num of levels in first group
    for(int first=k;first>=1;first--)
    {
        // cout << "\n\n\nfirst =" << first << endl;
        // cout << "k = " << k << endl;
        // cout << endl << endl;
        

        // levelsPerGroup = {a1, a2, ..., an} where ai = num of levels in group i
        vector<int> levelsPerGroup;

        levelsPerGroup.push_back(first); // Because first group has "first" many levels

        int count = first;  // Counts how many levels have been processed
        
        // For middle groups with k levels (it's possible there're no such groups)
        while(levels - count >= k)
        {
            count += k;
            levelsPerGroup.push_back(k);    
        }
        
        // For the last group with < k levels. if k levels, then handled above
        if(levels - count > 0)
            levelsPerGroup.push_back(levels - count);

        allGroupings.push_back(levelsPerGroup);

        // cout << "+++++++++++++++++++++++\n" << endl
        //      << "Size of First Window = "<< first << endl;
        // cout << "Set of Subset: " << Arrays::vecToStr(levelsPerGroup) << endl;

        if(k == levels) 
        {
            // cout << "\n\n\nwill break\n\n\n";
            break;
        }
    }

    return allGroupings;
}

vector<pair<int,int>> LevelGroupProcessor::prepareWithSubsets(int groupingNum, Config & c)
{
    vector<int> levelsPerGroup = allGroupings[groupingNum];
    vector<pair<int,int>> toFromPairs;
    
    Config allGroupWindows;
    
    int level = 1;
    for(int i=0; i<levelsPerGroup.size(); i++)
    {
        int from = level;
        level += levelsPerGroup[i];
        int to = level - 1;
        printf("Run Dynamic Programming from %d to %d \n", from, to);
        
        toFromPairs.push_back( pair<int,int>(to, from) );

        if(from == to)
            continue;;


        SingleRunConfig * data = chain->generateSubsetConfigFromBins(from, to);

        for(int x=0;x<10*10;x++)
            data->interLocs[0][x] = -9;

        int newLoc = 0;
        for(int x=0;x<data->numOfRNA;x++)
        {
            for(int y=x+1;y<data->numOfRNA;y++)
            {   
                int origID1 = data->rnaSequences[x]->origId;
                int origID2 = data->rnaSequences[y]->origId;
                int origLoc = weights->MatchingMatrixVal(origID1, origID2);
                if(origLoc == -9) 
                    continue;
                    // Because of even-even and odd-odd
                    //TODO: keep matrix even-odd only
                    
                data->rnaCollections[newLoc] = weights->GetWeightsTable(origLoc);
                data->interLocs[x][y] = newLoc;
                data->interLocs[y][x] = newLoc;
                // printf("x = %d, y = %d  --> newLoc = %d\n", x, y, newLoc);
                // printf("y = %d, x = %d  --> newLoc = %d\n", y, x, newLoc);
                newLoc++;
            
            }
        }
    
        // Do compression
        TempCompression tc;
        tc.doCompression(data, data->numOfRNA);

        PRBDPCore algo(data, 0, 1, loopAround);
        Config config = algo.config;
        // cout << "Config: " << endl;
        for(auto win : config)
        {
            // cout << win << "  " << win.weight << endl;
            allGroupWindows.add(win);
        }
        
        double e = config.totalWeight();
        cout << "Sub Energy of this Group = " << e << endl;

        delete data;
    }
    
    cout << "Total Energy of This Grouping = " << allGroupWindows.totalWeight() << endl;
    c = allGroupWindows;
    return toFromPairs;

}

Config LevelGroupProcessor::findMinEnergyStructure()
{
    double minEnergy = 0;
    Config minConfig;

    vector<vector<int> > x = GetAllGroupings();
    for(int groupNum=0; groupNum<x.size(); groupNum++)
    {
        Config c;
        prepareWithSubsets(groupNum, c);
        if(minEnergy > c.totalWeight())
        {
            minEnergy = c.totalWeight();
            minConfig = c;
        }
    }

    return minConfig;
}

