#include <iostream>
#include <cstdlib>
#include <ctime>
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

int main(int argc, char *argv[])
{
	srand(time(0));
	CmdLineArgs * args = Parser::ParseArgs(argc, argv);
	vector<OrigRNASeq> origRNASequences = Parser::GetAndParseInput();
	int numOfRNA = origRNASequences.size();

	Cache matchingCache;

	// string dir = "../rnaup_weights/output/";
	// string dir = "../output/";
	string dir = "output/";
	Weights weights(dir, &origRNASequences);
	weights.Read();

	// First find the initial, random structure
	// Chain chain = Chain::makeRandomChain(&origRNASequences);
	vector<int> permutation {3,1,0,2};
	Chain chain = Chain::makeGivenChain(&origRNASequences, permutation);
	chain.printFlatStruct();
	chain.determineStruct();

	if(args->k > chain.getBBHeight())
		args->k = chain.getBBHeight();
	
	LevelGroupProcessor lgProc(args->k, &chain, &weights, args->loopAround);
	Config minEnergyConfig = lgProc.findMinEnergyStructure();
	
	double lastRoundMinEnergy = minEnergyConfig.totalWeight();
	int rounds = 1;
	while(true)
	{
		cout << endl << "Round " << rounds << endl << endl;
		
		double thisRoundMinEnergy = 100000;
		Config thisRoundMinEnergyConfig;
		// Find nbr
		for(int start=0; start<numOfRNA-1; start++)
		{
			for(int end=start+1; end<numOfRNA; end++)	
			{
				// flip chain
				chain.flipSubChain(start, end);

				// cout << "start = " << start << "      end = " << end << endl << endl;

				LevelGroupProcessor lgProc(args->k, &chain, &weights, args->loopAround);
				Config config = lgProc.findMinEnergyStructure();				
				double energy = config.totalWeight();

				if(energy < thisRoundMinEnergy)
				{
					thisRoundMinEnergy = energy;
					thisRoundMinEnergyConfig = config;
				}

				
				// reset it to original state (at start of this round)
				chain.flipSubChain(start, end);
			}
		}

		// check if we made a better (lower) energy structure than before
		if(thisRoundMinEnergy < lastRoundMinEnergy)
		{
			// yes, we did, so continue finding better
			lastRoundMinEnergy = thisRoundMinEnergy;
			minEnergyConfig = thisRoundMinEnergyConfig;
		}
		else
		{
			break;
		}
	}

	cout << endl << endl;
	cout << "Last round min energy = " << lastRoundMinEnergy << endl;
	cout << "minEnergyConfig = " << minEnergyConfig << endl;



	return 0;
}



// vector<vector<int> > x = lgProc.GetAllGroupings();
	// for(int groupNum=0; groupNum<x.size(); groupNum++)
	// {
	// 	Config c;
 //        lgProc.prepareWithSubsets(groupNum, c);
 //        if(minEnergy > c.totalWeight())
 //        {
 //            minEnergy = c.totalWeight();
 //            minConfig = c;
 //        }
	// }