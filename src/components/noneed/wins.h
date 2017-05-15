typedef struct pair_pair 
{
    double G_is;
    int seq1_5,seq2_5,seq1_3,seq2_3;
    int id;
    int h,l;
} pair_pair;


typedef struct Event
{
    int t;
    int id;
    int type;
    pair_pair * pp;
} Event;


typedef struct argmin_obj 
{
    double val;
    int arg1;
    int arg2;
} argmin_obj;

