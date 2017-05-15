#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>
#include "Chain.h"
#include "rnaseq.h"
using namespace std;


unsigned long long rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned long long)hi << 32) | lo;
}



// Chain::Chain()
// {
//  head = createEmptyStructure();
// }

Chain::Chain(Pointer * h, vector<OrigRNASeq> * origRNASequences_) : head(h), origRNASequences(origRNASequences_)
{
    //throw std::runtime_error( "Constructor Chain(Pointer * h, vector<OrigRNASeq> * origRNASequences_) is not supported!" );
}   

Chain::Chain(vector<OrigRNASeq> * origRNASequences_, int numOfRNA) : 
            origRNASequences(origRNASequences_)
{
    head = makeSomeBinStruct(origRNASequences, numOfRNA);
}   

Chain::~Chain()
{
    destroyStruct(head);
}

Chain::Pointer * Chain::createEmptyStructure()
{
    Chain::Pointer * head = new Chain::Pointer;
    head->elem = NULL;
    return head;
}

Chain::Pointer * Chain::makeSomeBinStruct(vector<OrigRNASeq> * origRNASequences, int numOfRNA)
{
    Chain::Pointer * head = createEmptyStructure();
    
    Chain::EOB * n = new Chain::EOB;
    n->next = NULL;
    n->origId = (origRNASequences->at(0)).origId;
    n->type = (origRNASequences->at(0)).type;
        
    head->elem = n;
    
    for(int i=1;i<numOfRNA;i++)
    {
        Chain::EOB * n2 = new Chain::EOB;
        n2->next = NULL;
        n2->origId = (origRNASequences->at(i)).origId;
        n2->type = (origRNASequences->at(i)).type;
        n->next = n2;       
        n = n2;
    }
    
    return head;
}

Chain Chain::makeRandomChain(vector<OrigRNASeq> * origRNASequences)
{
    vector<int> v;
    for(int i=0; i< origRNASequences->size(); i++)
        v.push_back(i);

    // shuffle
    int N = v.size();
    for(int i=0; i<N; i++)
    {
        int randPos = rand() % N;
        int temp = v[i];
        v[i] = v[randPos];
        v[randPos] = temp;
    }

    // cout << "Random chain:" << endl;
    // for(auto e : v)
    //     cout << e << " ";
    // cout << endl;

    return Chain::makeGivenChain(origRNASequences, v);
}

Chain Chain::makeGivenChain(vector<OrigRNASeq> * origRNASequences, vector<int> sequence)
{
    // A quick hack to prevent this function from going down
    // if an empty vector is passed -- shouldn't happen ideally
    if(sequence.size() == 0)
        sequence.push_back(0);

    Chain::Pointer * head = createEmptyStructure();
    
    int pos = sequence[0];
    
    Chain::EOB * n = new Chain::EOB;
    n->next = NULL;
    n->origId = (origRNASequences->at(pos)).origId;
    n->type = (origRNASequences->at(pos)).type;
        
    head->elem = n;
    
    for(int i=1;i<sequence.size();i++)
    {
        pos = sequence[i];
        Chain::EOB * n2 = new Chain::EOB;
        n2->next = NULL;
        n2->origId = (origRNASequences->at(pos)).origId;
        n2->type = (origRNASequences->at(pos)).type;
        
        n->next = n2;
        
        n = n2;
    }
    
    return Chain(head, origRNASequences);
}

void Chain::destroyStruct(Pointer * head)
{
    EOB * n = head->elem;
    delete head;
    destroyStruct(n, 0);
}

void Chain::destroyStruct(EOB * n, int i)
{
    if(n->next == NULL)
    {
        delete n;
        return;
    }
    
    destroyStruct(n->next, i+1);
    delete n;
}

void Chain::flipSubChain(int start, int end)
{
    EOB * n = head->elem;
    EOB * startEOB, * endEOB;
    
    int i = 0;
    EOB * array[20];
    EOB * prevToStart = NULL;
    int flipCount = 0;
    while(n != NULL)
    {   
        if(i == start - 1)
            prevToStart = n;
            
        if(i == start)
        {
            startEOB = n;
            // printf("Start id = %d \n", startEOB->origId);
        }   
        else if(i == end)
        {
            endEOB = n;
            // printf("End id = %d \n", endEOB->origId);
        }
        
        if(i >= start && i <= end)
        {
            //printf("reversing pointers of %d \n", n->origId);
            
            array[flipCount++] = n;
        }
            
        n = n->next;
        i++;
    }
    
    EOB * nextToEnd = endEOB->next;
    
    for(i=flipCount-1;i>0;i--)
    {
        EOB * thisOne = array[i];
        EOB * realPrev = array[i-1];
        
        thisOne->next = realPrev;
    }
    
    startEOB->next = nextToEnd;
    if(prevToStart == NULL) //, then prev is head pointer;
    {
        head->elem = endEOB;
    }
    else
        prevToStart->next = endEOB;
    
    
    // printf("reversed... \n");
}

void Chain::printChain()
{

    EOB * n = head->elem;
    int c = 0;
    printf("Chain: \n");
    while(n != NULL)
    {   
        if(n->type == 2)
            printf("| ");
        else if(n->type == 0)
            printf("e%d ", n->origId);
        else if(n->type == 1)
            printf("o%d ", n->origId);
            
        n = n->next;
        c++;
    }
    printf("\n\n");
}

void Chain::determineStruct()
{
    EOB * n = head->elem;
    
    int currKnown = n->type;
    int level = 1;
    printf("Chain: \n");
    printf("Level 1: ");
    while(n != NULL)
    {   
        currKnown = n->type;
    
        if(n->type == 2)
            printf("| ");
        else if(n->type == 0)
            printf("e%d ", n->origId);
        else if(n->type == 1)
            printf("o%d ", n->origId);
            
        n = n->next;
        
        if(n != NULL)
            if(currKnown != n->type)
            {
                printf("\nLevel %d: ", ++level);
                //printf("\n");
            }
    }
    printf("\n\n");

}

int Chain::getBBHeight()
{
    EOB * n = head->elem;
    int currKnown = n->type;
    int height = 1;
    while(n != NULL)
    {   
        currKnown = n->type;
        n = n->next;
        
        if(n != NULL)
            if(currKnown != n->type)
            {
                height++;
            }
    }
    return height;
}

SingleRunConfig * Chain::generateSubsetConfigFromBins(int from, int to)
{
    EOB * n0 = head->elem;

    int currKnown = n0->type;
    int level = 1;
    EOB * fromNode = NULL;
    while(n0 != NULL)
    {   

        currKnown = n0->type;

        if(level == from)
        {
            fromNode = n0;
            break;
        }

        n0 = n0->next;

        if(n0 != NULL)
            if(currKnown != n0->type)
                ++level;

    }
    
    SingleRunConfig * data = new SingleRunConfig;
    
    for(int i=0;i<10;i++) 
        data->newId2OldId[i] = -9;
    for(int i=0;i<10*10;i++) 
    {
        data->setup[0][i] = -9;
        data->parentOf[0][i] = -9;
    }
    
    //to create setup[][] array, just go thorugh the sequences linearly
    
    //as in for(i=from;i<=to;i++)
    EOB * n = n0, * firstNodeLastLevel = NULL;
    int num = 0, setupRow = 0, setupColumn = 0;
    
    string stringRep;
    
    while(level <= to && n != NULL)
    {
        currKnown = n->type;
        
        data->newId2OldId[num] = n->origId;
        data->setup[setupRow][setupColumn] = num;

        stringRep += std::to_string(n->origId);
        
        RNASeq * rnaSeq = new RNASeq;
        rnaSeq->origId = n->origId;
        rnaSeq->id = num;
        rnaSeq->type = (origRNASequences->at(n->origId)).type;
        rnaSeq->string = (origRNASequences->at(n->origId)).string;
        rnaSeq->originalLength = (origRNASequences->at(n->origId)).originalLength;
        rnaSeq->name = (origRNASequences->at(n->origId)).name;
        data->rnaSequences.push_back(rnaSeq);

        num++;

        if(level == to and firstNodeLastLevel == NULL)
            firstNodeLastLevel = n;
                    
        n = n->next;

        setupColumn++;

        // The following block detects a new level in the chain
        // and changes the counters to that of a new level
        if(n != NULL)
        {
            if(currKnown != n->type)
            {
                ++level;
                setupRow++;
                setupColumn = 0;
                
                if(level <= to)
                    stringRep += "|";
            }
            else
                stringRep += "_";
        }
    }

    data->numOfRNA = num;   
    data->totalLevels = to - from + 1;
    data->stringRep = stringRep;

    //fix parentOf matrix
    for(int i=0;i<9;i++)    //rows
    {
        for(int j=0;j<10;j++)   //cols on top level
        {
            if(data->setup[i][j] == -9) continue;
            
            int parentRNA = data->setup[i][j];
            for(int k=0;k<10;k++)   //cols on second level
            {
                if(data->setup[i+1][k] == -9) continue;
                
                int childRNA = data->setup[i+1][k];
                data->parentOf[parentRNA][childRNA] = 1;
                
                // printf("setting parent=1 for %d > %d \n", parentRNA, childRNA);
            }
        }
    }

    // Wraparound: Last row elems are parent of the first row
    // However, we can wrap around only if first,last are even,odd or odd,even
    // i.e, the number of levels is even
    if(data->totalLevels >= 4 and data->totalLevels % 2 == 0)
    {
        int i = data->totalLevels-1;
        printf("i = %d\n", i);
        for(int j=0;j<10;j++)   //cols on top level
        {
            if(data->setup[i][j] == -9) continue;   // if no RNA here, do nothing
            
            int parentRNA = data->setup[i][j];
            for(int k=0;k<10;k++)   //cols on second level
            {
                if(data->setup[0][k] == -9) // Note: 0 for first level
                    continue;
                
                int childRNA = data->setup[0][k]; // Note: 0 for first level
                data->parentOf[parentRNA][childRNA] = 1;
                
                // printf("setting parent=1 for %d > %d \n", parentRNA, childRNA);
            }
        }
    }


    // printf("exiting...\n");

    
    return data;
}

void Chain::printFlatStruct2File(std::FILE * file)
{
    EOB * n = head->elem;
    
    int currKnown = n->type;
    // int level = 1;
    while(n != NULL)
    {   
        currKnown = n->type;
    
        if(n->type == 0)
            fprintf(file, "e%d ", n->origId);
        else if(n->type == 1)
            fprintf(file, "o%d ", n->origId);
            
        n = n->next;
        
        if(n != NULL)
            if(currKnown != n->type)
            {
                fprintf(file, "| ");
            }
    }
}

void Chain::printFlatStruct()
{
    EOB * n = head->elem;
    
    int currKnown = n->type;
    // int level = 1;
    while(n != NULL)
    {   
        currKnown = n->type;
    
        if(n->type == 0)
            printf("e%d ", n->origId);
        else if(n->type == 1)
            printf("o%d ", n->origId);
            
        n = n->next;
        
        if(n != NULL)
            if(currKnown != n->type)
            {
                printf("| ");
            }
    }
    std::cout << std::endl;
}




// Pointer * Chain::makeRandomBinStruct(OrigRNASeq ** origRNASequences, int numOfRNA)
// Chain Chain::makeRandomBinStruct(OrigRNASeq ** origRNASequences, int numOfRNA)
// {

    // //srand(time(NULL));
    // srand(rdtsc());
    // int array[MAX_NUM_RNA];
    // int i;
    
    // for(i=0;i<numOfRNA;i++)
    // {
    //  array[i] = origRNASequences[i]->origId;
    // }
    // //now randomize it
    // for(i=0;i<numOfRNA;i++)
    // {
    //  //select index 2
    //  int r1 = rand() % numOfRNA; //0,1,2,..,n-1

    //  int temp = array[i];
    //  array[i] = array[r1];
    //  array[r1] = temp;
    // }
    
    
    // Pointer * head = createEmptyStructure();
    
    // EOB * n = (EOB*) malloc(sizeof(EOB));
    // n->next = NULL;
    // n->origId = origRNASequences[array[0]]->origId;
    // //printf("n->origId = %d \n", n->origId);
    // n->type = origRNASequences[array[0]]->type;
        
    // head->elem = n;
    
    // //printf("id: %d \n", head->elem->origId);
    // for(i=1;i<numOfRNA;i++)
    // {
    //  EOB * n2 = (EOB*) malloc(sizeof(EOB));
    //  n2->next = NULL;
    //  n2->origId = origRNASequences[array[i]]->origId;
    //  //printf("n2->origId = %d \n", n2->origId);
    //  n2->type = origRNASequences[array[i]]->type;
        
    //  //printf("id: %d \n", n2->origId);
        
    //  n->next = n2;
        
    //  n = n2;
    // }
    
    // return Chain(head,orig);
    // return Chain();
// }


// int main()
// {
//  Chain c;

//  return 0;
// }