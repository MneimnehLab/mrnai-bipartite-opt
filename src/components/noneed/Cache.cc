#include <vector>
#include "Cache.h"
using namespace std;


Cache::Cache()
{

}

// void Cache::AddSubEnsemble(SubEnsemble * subensemble)
void Cache::AddConfig(string stringRep, Config config)
{
    // matching_hashtable[stringRep] = subensemble;
    matching_hashtable[stringRep] = config;
}

// SubEnsemble * Cache::GetEnsemble(string stringRep)
Config Cache::GetConfig(string stringRep)
{
    // unordered_map<string, SubEnsemble*>::const_iterator elemLoc 
    //  = matching_hashtable.find (stringRep);
    // if ( elemLoc == matching_hashtable.end() )
    //  return nullptr;

    // return (*elemLoc).second;

    unordered_map<string, Config>::const_iterator elemLoc 
        = matching_hashtable.find (stringRep);
    if ( elemLoc == matching_hashtable.end() )
        return Config();

    return (*elemLoc).second;   
}



    