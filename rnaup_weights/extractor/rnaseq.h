#define MAX_RNA_SIZE 250
#define MAX_NUM_RNA 10
#define GAP 4

typedef struct OrigRNASeq
{
	int origId;
	char * string;
	int type;
	int originalLength;
	char name[6];
	
} OrigRNASeq;

typedef struct RNASeq
{
	int id;
	char * string;
	char name[6];
	
	int * compressedRNA;
	int * expandedRNAmap;
	int interactsWith[10];
	
	//int * parentOf[10];
	
	int compressedLength;
	int originalLength;
	
	int type;	//even or odd
	
	int origId;
	
} RNASeq;



typedef struct InteractWeight
{
	double **** weights;
} InteractWeight;



typedef struct SingleRunConfig
{
	int numOfRNA;
	char ** strArray1;
	int setup[10][10];
	int parentOf[10][10];
	int interLocs[10][10];
	int totalLevels;
	double **** rnaCollections[50];
	int rnaCollSize;
	int newId2OldId[10];
	char * stringRep;
	int gapSize;

	RNASeq * rnaSequences[10];
} SingleRunConfig;

typedef struct InterDims
{
	int n1;
	int n2;
} InterDims;



typedef struct CmdLineArgs
{
	int  parallel;
	int  k;
	int  fillGaps;
	int  trials;
	int  gapSize;
	int  winSize;
	int  GU;
	char rnaupOut[50];
	char fileName[50];
} CmdLineArgs;



typedef struct Matching
{
	double energy;
	
	
} Matching;



/*
typedef struct OneRNAup
{
	double
} OneRNAup;
*/

typedef struct OneInteractionStructure
{
	int topSeq, botSeq, startOnTopSeq, indexOnTopSeq, startOnBotSeq, indexOnBotSeq;
	int trueTopSeq, trueBotSeq;
	double weight;
} OneInteractionStructure;

typedef struct SubEnsemble
{
	OneInteractionStructure * interactions[100];
	int count;
	double subenergy;
	char * stringRep;
	int firstRowRNAs[MAX_NUM_RNA], lastRowRNAs[MAX_NUM_RNA];
	
} SubEnsemble;

typedef struct Ensemble
{
	int k;
	int first;
	int count;
	double energy;
	SubEnsemble * subs[10];
} Ensemble;



typedef struct RNAupOutput
{
	
} RNAupOutput;
