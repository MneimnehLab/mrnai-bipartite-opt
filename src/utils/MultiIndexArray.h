#ifndef _MULTI_INDEX_ARRAY_H
#define _MULTI_INDEX_ARRAY_H

#include <vector>
// #include <iostream>
#include <initializer_list>


// e.g. array3D[0][0][2] <=> container[2]
template<typename T>
class MultiIndexArray
{

private:
    // T * container;
    std::vector<T> container;
    std::vector<int> dimSizes;
    std::vector<int> bases;
    int totalDimSize;
    int totalLevels;

    inline int doHash(const int * arr);

public:
    //@param dimSizes A list of the element in each dimension
    // Example [4, 5, 2] is a 4x5x2  3-D array 
    // MultiIndexArray(std::initializer_list<int> l);
    MultiIndexArray();
    MultiIndexArray(std::vector<int> v);
    // ~MultiIndexArray();
    // MultiIndexArray(MultiIndexArray<T>&& other);
    // MultiIndexArray<T>& operator=(MultiIndexArray<T>&& other);

    T& operator[] (const int* index);
    T& operator[] (const int index);

    // static MultiIndexArray CreateMultiIndexArray(std::initializer_list<T> l);
    static MultiIndexArray  CreateMultiIndexArray(std::vector<int> v);

    inline void revHash(int hash, int * x);
    inline vector<int> revHash(int hash);

};


template<typename T>
MultiIndexArray<T>::MultiIndexArray()
{
    totalLevels = 1;
    totalDimSize = 1;
    
    // make hash multipliers, i.e., the bases of the hash function. bases[0] = len(Level1), base[1] = base[0]*len(level2), and so on...
    // (ADD more details in hash function comments...)
    bases.resize(1);
    dimSizes.resize(1);
    bases[0] = dimSizes[0] = 1;
    
    container.resize(1);
}


template<typename T>
// MultiIndexArray<T>::MultiIndexArray(std::initializer_list<int> l): dimSizes(l)
MultiIndexArray<T>::MultiIndexArray(std::vector<int> l): dimSizes(l)
{
    totalLevels = dimSizes.size();
    totalDimSize = 1;
    for(int i=0;i<totalLevels;i++)
    {
        totalDimSize *= dimSizes[i];
        // std::cout << dimSizes[i] << std::endl;
    }
    
    // make hash multipliers, i.e., the bases of the hash function. bases[0] = len(Level1), base[1] = base[0]*len(level2), and so on...
    // (ADD more details in hash function comments...)
    bases.resize(totalLevels);
    bases[totalLevels-1] = dimSizes[totalLevels-1];
    for(int i=totalLevels-2;i>=0;i--)
        bases[i] = dimSizes[i] * bases[i+1];

    // container = new T[totalDimSize];
    // container.resize(totalDimSize);
    container = std::vector<T>(totalDimSize, T());
}


template<typename T>
int MultiIndexArray<T>::doHash(const int * arr)
{
    int sum = arr[totalLevels-1];
    for(int i=totalLevels-2;i>=0;i--)
        sum += arr[i] * bases[i+1];
    return sum;
}


template<typename T>
void MultiIndexArray<T>::revHash(int hash, int * x)
{
    x[totalLevels-1] = hash % dimSizes[totalLevels-1];
    for(int i=totalLevels-2;i>=0;i--)
    {
        x[i] = (hash / bases[i+1]) % dimSizes[i];
    }
}

template<typename T>
vector<int> MultiIndexArray<T>::revHash(int hash)
{
    vector<int> v(totalLevels);
    v[totalLevels-1] = ( hash % dimSizes[totalLevels-1] ) ;
    for(int i=1;i<totalLevels;i++)
    {
        v[i] = (  (hash / bases[i+1]) % dimSizes[i] );
    }
    return v;
}

template<typename T>
T& MultiIndexArray<T>::operator[] (const int* indices)
{
    int index = doHash(indices);
    return container[index];
}

template<typename T>
T& MultiIndexArray<T>::operator[] (const int index)
{
    return container[index];
}


// template<typename T>
// MultiIndexArray MultiIndexArray<T>::CreateMultiIndexArray(std::initializer_list<T> l)
// {

// }
    
template<typename T>
MultiIndexArray<T>  MultiIndexArray<T>::CreateMultiIndexArray(std::vector<int> v)
{
    return MultiIndexArray<T>(v);
}

// template<typename T>
// MultiIndexArray<T>::~MultiIndexArray()
// {
//  cout << "MultiIndexArray Destructor called" << endl;
//  delete container;
// }

// template<typename T>
// MultiIndexArray<T>::MultiIndexArray(MultiIndexArray<T>&& other)
// {
//  cout << "Move constructor called" << endl;
//  container = other.container;
//  other.container = nullptr;
//  dimSizes = other.dimSizes;
//  bases = other.bases;
//  totalDimSize = other.totalDimSize;
//  totalLevels = other.totalLevels;

// }

// template<typename T>
// MultiIndexArray<T>& MultiIndexArray<T>::operator=(MultiIndexArray<T>&& other)
// {
//  cout << "Move assignment operator called" << endl;
//  container = other.container;
//  other.container = nullptr;
//  dimSizes = other.dimSizes;
//  bases = other.bases;
//  totalDimSize = other.totalDimSize;
//  totalLevels = other.totalLevels;

// }

#endif
