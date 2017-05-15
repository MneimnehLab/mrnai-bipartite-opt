#include <iostream>
#include <sstream>
#include <fstream>
#include "Weights.h"
#include "RNAProperties.h"
#include "../utils/Arrays.h"
using std::string;

RNAProperties::RNAProperties(Weights * weights, vector<OrigRNASeq> * v) :
	weights(weights), origRNASequences(v)
{
	// separate RNAs into even and odd vectors,
    // depending on labelling
    for(auto seq : *origRNASequences)
    {
        if(seq.type) // odd (== 1)
            oddSeq.push_back(seq.origId);
        else
            evenSeq.push_back(seq.origId);

        rnaLengths.push_back(seq.originalLength); 
    }

    numEven = evenSeq.size();
    numOdd = oddSeq.size();
}


int RNAProperties::getNumEven()
{
    return numEven;
}

int RNAProperties::getNumOdd()
{
    return numOdd;
}

Weights * RNAProperties::getWeights()
{
	return weights;
}
