#ifndef _CACHE_H
#define _CACHE_H

#include <unordered_map>
#include <string>
#include "rnaseq.h"
#include "Config.h"
using namespace std;


class Cache
{
private:
    // unordered_map<string, SubEnsemble*> matching_hashtable;
    unordered_map<string, Config> matching_hashtable;

public:
    Cache();
    // void AddSubEnsemble(SubEnsemble * subensemble);
    // SubEnsemble * GetEnsemble(string stringRep);
    void AddConfig(string, Config);
    Config GetConfig(string);
};

#endif
