//Even_Odd_Barrier :- EOB
//typedef struct EOB EOB;
//typedef struct Pointer Pointer;

typedef struct EOB
{
	int type;	//0: even, 1: odd
	int origId;
	
	struct EOB * next;
} EOB;

typedef struct Pointer
{
	EOB * elem;
} Pointer;
