	/*
  Last changed Time-stamp: <2008-07-04 16:15:50 ulim>
  $Id: RNAup.c,v 1.5 2008/07/04 14:27:09 ivo Exp $
  
  Ineractive Access to cofolding routines
  c Ivo L Hofacker
  Vienna RNA package
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <float.h> 
#include <time.h> 
#include <stdarg.h>
#include "fold.h"
#include "fold_vars.h"
#include "PS_dot.h"
#include "utils.h"
#include "part_func.h"
#include "part_func_up.h"
#include "duplex.h"
#include "energy_const.h"
#include "../extractor/wins.h"
#include "../extractor/sched.c"
#include "../extractor/rnaseq.h"

extern long startT;
extern long lastT;


extern void read_parameter_file(const char fname[]);


/*@unused@*/
static char rcsid[] = "$Id: RNAup.c,v 1.5 2008/07/04 14:27:09 ivo Exp $";

#define PRIVATE static
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define EQUAL(A,B) (fabs((A)-(B)) < 1000*DBL_EPSILON)
PRIVATE void tokenize(char *line, int * numOfRNA, char **strArray1);
PRIVATE char *tokenize_one(char *line);
PRIVATE int comp_nums(const int *num1, const int *num2);
PRIVATE int get_u_values(char unstrs[], int **u_vals, int l1);
PRIVATE void seperate_bp(char **inter, int len1, char **intra_l, char **intra_s);
PRIVATE void print_interaction(interact *Int, char *s1, char *s2, pu_contrib *p_c, pu_contrib *p_c2, int w, int incr3, int incr5);
PRIVATE void print_unstru(pu_contrib *p_c, int w);

//int doWins(char* str1, char* str2);
void showChildrenForW(int w, int count, FILE*);
int makeListForWk(char* str1, char* str2, int k);
void filterAllW(char* string1, char* string2);

void perform2(char * string1, char * string2, int ic, double ***** rnaupCollections);
void reverseRNAupOutput(int location, int n1, int n2, double ***** rnaupCollections);	//rnaWins

static char scale1[] = "....,....1....,....2....,....3....,....4";
static char scale2[] = "....,....5....,....6....,....7....,....8";

PRIVATE void usage(void);
/* defaults for -u and -w */
PRIVATE int default_u; /* -u options for plotting: plot pr_unpaired for 4 nucleotides */
PRIVATE int default_w; /* -w option for interaction: maximal region of interaction is 25 nucleotides */
PRIVATE double RT;
extern double INF1;

//extern double G_holder[120][120][120][120];
//extern double rnaWins[120][120][26];
extern double **** rnaWins;
pair_pair ** ppList;


char *dummy = NULL, *temp = NULL, *line = NULL;
char *structure = NULL, *cstruc = NULL, *cstruc_l = NULL, *cstruc_s = NULL;
char fname[53], ffname[53], temp_name[201], first_name[53], my_contrib[10];
char up_out[250], unstrs[201], name[400], cmd_line[500];
char *ParamFile = NULL;
char *ns_bases = NULL, *c, *head;
int i, length1, length2, length, l, sym, r, Switch, header, output;
double energy, min_en;
double sfact = 1.07;
int istty;
int noconv = 0;
/* variables for output */
pu_contrib *unstr_out, *unstr_short;
interact *inter_out;
/* pu_out *longer; */
char *title;
/* commandline parameters */
int w; /* length of region of interaction */
int incr3; /* add x unpaired bases after 3'end of short RNA*/
int incr5; /* add x unpaired bases after 5'end of short RNA*/
int unstr; /* length of unpaired region for output*/
int upmode; /* 1 compute only pf_unpaired, >1 compute interactions 
	  2 compute intra-molecular structure only for long RNA, 3 both RNAs */
int task; /* input mode for calculation of interaction */


double matrix[MAX_RNA_SIZE][MAX_RNA_SIZE][26][26];
char stringEven[MAX_RNA_SIZE];
char stringOdd[MAX_RNA_SIZE];
InterDims interDimArray[10];

char * rnaupOut;

void initRNAup()
{
	noGU = 1;
	printf("Our Version!\n");
    

	RT = ((temperature + K0) * GASCONST / 1000.0);
}

void doRNAupAndReverse(char * evenString, char * oddString, int location, double ***** rnaupCollections, char * rnaupOut1, int GU)
{
	// printf("evenString = %s\n", evenString);
	noGU = 1 - GU;
	//char * ns_bases = "AC,GA,CA,UG";
	char * ns_bases = NULL;
	// Non canonical pairs 
	if (ns_bases != NULL)
    {
        nonstandards = space(33);
        c = ns_bases;
        i = sym = 0;
        if (*c == '-')
        {
            sym = 1;
            c++;
        }
        while (*c != '\0')
        {
            if (*c != ',')
            {
                nonstandards[i++] = *c++;
                nonstandards[i++] = *c;
                if ((sym) && (*c != *(c - 1)))
                {
                    nonstandards[i++] = *c;
                    nonstandards[i++] = *(c - 1);
                }
            }
            c++;

        }
    }
	
	// End Non canonical pairs
	
	rnaupOut = rnaupOut1;

	// printf("Starting RNAup... ");
	// printMilliSecs();

	strcpy(stringEven, evenString);
	strcpy(stringOdd, oddString);
			
	// printf("stringEven = %s \n", stringEven);
	// printf("stringOdd = %s \n", stringOdd);
			
	int n1, n2;
	n1 = strlen(stringEven);
	n2 = strlen(stringOdd);
			
			
	// printf("Sending to perform2... ");
	// printMilliSecs();

	perform2(stringEven, stringOdd, location, rnaupCollections);


	// printf("Exited perform2... ");
	// printMilliSecs();

	
	// printf("Sending to reverse... ");
	// printMilliSecs();

	reverseRNAupOutput(location, n1, n2, rnaupCollections);

	// printf("Got back from reverse... ");
	// printMilliSecs();
/*	
	int x,y,w;
	for(x=0;x<50;x++)
		for(y=0;y<50;y++)
			for(w=1;w<=25;w++)
				printf("rnaWins[%d][%d][%d] = %f \n", x,y,w, rnaWins[x][y][w]);
*/	


}





void perform2(char * string1, char * string2, int ic, double ***** rnaupCollections)
{
//	int *** returnPointer;

	RT = ((temperature + K0) * GASCONST / 1000.0);
	cut_point=-1;
	//printf("\n\n--------------------------------------------------\nperform2 \n");
	/* default settings for RNAup */
    head = NULL; /* header text - if header wanted, see header */
    header = 1; /* if header is 0 print no header in output file: option -nh */
    output = 1; /* if output is 0 make no output file: option -o */
    Switch = 0; /* the longer sequence is selected as the target */
    task = 0;
    upmode = 1; /* default is one sequence, option -X[p|f] has to be set
		 for the calculation of an interaction, if no "&" is in
		 the sequence string  */
    unstrs[0] = '\0';
    default_u = 4;
    unstr = default_u;
    default_w = 25;
    w = default_w;
    
    incr3 = 0;
    incr5 = 0;
    do_backtrack = 1;
    length1 = length2 = 0;
    title = NULL;
    unstr_out = NULL;
    inter_out = NULL;
    my_contrib[0] = 'S';
    my_contrib[1] = '\0';
    first_name[0] = '\0';
	
	
		int * u_vals = NULL;
		
		//printf("in RNAup, n1 = %d, n2 = %d \n", strlen(string1), strlen(string2));
		//printf("RNA1: %s \n", string1);
		//printf("RNA2: %s \n", string2);
		
		/* main loop: continue until end of file */
		//printf("count = %d \n", ic);
        if (string1 != NULL)
        {
            length1 = (int) strlen(string1);
        } 
		else
        {
            nrerror("sequence is NULL, check your input.");
        }
        
		if (string2 != NULL)
		{
			length2 = (int) strlen(string2);
		}
		else
		{
			nrerror("one of the sequences is NULL, check your input.");
		}

		/* write longer seq in string1 and and shorter one in string2 */
		if (length1 < length2 && Switch)
		{
			length = length1;
			length1 = length2;
			length2 = length;

			temp = (char *) space(sizeof (char) *strlen(string1) + 1);
			(void) sscanf(string1, "%s", temp);
			string1 = (char *) xrealloc(string1, sizeof (char) *length1 + 1);
			(void) sscanf(string2, "%s", string1);
			string2 = (char *) xrealloc(string2, sizeof (char) *length2 + 1);
			(void) sscanf(temp, "%s", string2);

			free(temp);
			temp = NULL;
			
		}
	
	
		/* get values for -u */
        if (!get_u_values(unstrs, &u_vals, length1))
        {
            nrerror("option -u: length value exceeds sequence length\n");
        }

        for (l = 0; l < length1; l++)
        {
            string1[l] = toupper(string1[l]);
            if (!noconv && string1[l] == 'T') string1[l] = 'U';
        }
        for (l = 0; l < length2; l++)
        {
            string2[l] = toupper(string2[l]);
            if (!noconv && string2[l] == 'T') string2[l] = 'U';
        }
        if (length1 > length2)
        {
            structure = (char *) space(sizeof (char) *(length1 + 1));
        } 
		else
        {
            structure = (char *) space(sizeof (char) *(length2 + 1));
        }
        update_fold_params();
        if (cstruc_s != NULL)
            strncpy(structure, cstruc_s, length2 + 1);
		/*printf("fold 1 \n");
		printf("string1: %s\n", string1);
		printf("structure: %s\n", structure);*/
        min_en = fold(string1, structure);
        (void) fflush(stdout);
		int wplus, w_sh;
		
		/* calculate prob. unstruct. for shorter seq */
		w_sh = w;
		/* len of unstructured region has to be <= len shorter seq. */
		if (w > length2) w_sh = length2;
		if (cstruc_s != NULL)
			strncpy(structure, cstruc_s, length2 + 1);
		//printf("fold 2 \n");
		min_en = fold(string2, structure);
		pf_scale = exp(-(sfact * min_en) / RT / length2);
		
		
		if (length2 > 2000) fprintf(stderr, "scaling factor %f\n", pf_scale);
		init_pf_fold(length2);
		if (cstruc_s != NULL)
			strncpy(structure, cstruc_s, length2 + 1);
		
		
		energy = pf_fold(string2, structure);
		
		/*
		
		char sequenceXYZ[] = "ACAAUACAGAGAUGAUCAGCAGUUCCCCUGCAUAAGGAUGAACCGUUUUACAAAGACCAUUGCACUCCGGUCUUUGCCUUUUGGCUUAGAUCAAGUGUAGUAUCUGU";
		char * structureXYZ = (char*)malloc(sizeof(char) * strlen(sequenceXYZ));
		//energy = pf_fold(string2, structure);
		energy = pf_fold(sequenceXYZ, structureXYZ);
		printf("\n\npf_fold for string2 = %f\n", energy);
		printf("structure... = %s\n\n", structureXYZ);
		
		exit(0);
		
		*/
		
		unstr_short = pf_unstru(string2, w_sh);
		free_pf_arrays(); /* for arrays for pf_fold(...) */
	
		/* calculate prob. unstructured for longer seq */
		wplus = w + incr3 + incr5;
		/* calculate prob. unpaired for the maximal length of -u */
		if (u_vals[u_vals[0]] > wplus) wplus = u_vals[u_vals[0]];
		/* length of the unstructured region has to be <= len longer seq. */
		if (wplus > length1) wplus = length1;
		if (cstruc_l != NULL)
			strncpy(structure, cstruc_l, length1 + 1);
			

			//printf("fold 3 \n");
		min_en = fold(string1, structure);
		pf_scale = exp(-(sfact * min_en) / RT / length1);
		
		
		
		if (length1 > 2000) fprintf(stderr, "scaling factor %f\n", pf_scale);
		init_pf_fold(length1);
		if (cstruc_l != NULL)
			strncpy(structure, cstruc_l, length1 + 1);
		
		energy = pf_fold(string1, structure);
		unstr_out = pf_unstru(string1, wplus);
		
		/** At this point, we have P_u[i,j] fora ll i,j s.t. j-i < 25 */
					
		free_pf_arrays(); /* for arrays for pf_fold(...) */
		/* now make output to stdout and to the output file */

		/* calculate interaction between two sequences */
		int count;



		//ADD STACK call here


		//printf("about to send to pf_interact for n1=%d, n2=%d \n", strlen(string1), strlen(string2));	
		inter_out = pf_interact(string1, string2, unstr_out, unstr_short, w, cstruc, incr3, incr5);
		print_interaction(inter_out, string1, string2, unstr_out, unstr_short, w, incr3, incr5);
		

		int n1a = strlen(string1);
		int n2a = strlen(string2);

		rnaupCollections[ic] = rnaWins;
		

		
/*		printf("(*rnaupCollections)[1][1][1] = %f \n",(*rnaupCollections)[1][1][1]);
exit(0);
*/
/*
	int x,y,w;
	for(x=1;x<50;x++)
		for(y=1;y<50;y++)
			for(w=1;w<=25;w++)
				printf("rnaupCollections[ic][%d][%d][%d] = %f \n", x,y,w, rnaupCollections[ic][x][y][w]);
//				printf("rnawins[%d][%d][%d] = %f \n", x,y,w, rnaWins[x][y][w]);
*/

		
		//printf("here 3 \n");
        if (structure != NULL) free(structure);
        structure = NULL;
//printf("here 4 \n");
        if (title != NULL) free(title);
        title = NULL;
//printf("here 5 \n");
        if (u_vals != NULL) free(u_vals);
        u_vals = NULL;
//printf("here 6 \n");
		free_pu_contrib(unstr_out);
		free_interact(inter_out);
        free_pu_contrib(unstr_short);
        free_arrays(); /* for arrays for fold(...) */
//printf("here 7 \n");
        if (cstruc != NULL) free(cstruc);
        cstruc = NULL;

        if (cstruc_l != NULL) free(cstruc_l);
        cstruc_l = NULL;

        if (cstruc_s != NULL) free(cstruc_s);
        cstruc_s = NULL;

        (void) fflush(stdout);
//printf("here 8 \n");
        /*if (string1 != NULL)
        {
            free(string1);
            string1 = NULL;
        }
        if (string2 != NULL) free(string2);
        string2 = NULL;*/
//printf("here 9 \n");		
		if (line != NULL) free(line);
		//if (string1 != NULL) free(string1);
		//if (string2 != NULL) free(string2);
		if (cstruc != NULL) free(cstruc);
		if (cstruc_l != NULL) free(cstruc_l);
		if (cstruc_s != NULL) free(cstruc_s);
//printf("here 10 \n");	

//	return 0;
	
	
}









































PRIVATE void usage(void)
{
    nrerror("usage:\n"
            "RNAup [-u list] [-w len] [-b] [-Xp|-Xf] [-c \"SHIME\"] [-5 incr]\n"
            "      [-3 incr] [-target] [-o] [-C] [-T temp] [-noLP]\n"
            "      [-d[0|2]] [-noGU] [-noCloseGU] [-P paramfile] [-4]\n"
            "      [-nsp pairs] [-S scale] [-noconv] \n");
}

/* call:  tokenize(line,&seq1,&seq2); the sequence string is split at the "&"
   and the first seq is written in seq1, the second into seq2  */



   
   
   
/* using sscanf instead of strcpy get's rid of trainling junk on the input line */
void tokenize(char *line, int * numOfRNA, char **strArray1)
{
    //char * line = "RSTUVWXYZ&LMNO&12345&ABCDEFGHIJK";
	
	char *start = line, *pos = line;
    int cut = -1;
    int i = 0 ;
	
	while(pos = strchr(pos + 1, '&'))
	{
		cut = (int) (pos - start) + 1;
		char * seq1 = (char *) space((cut + 1) * sizeof (char));	
		strncpy(seq1, start, cut);
		seq1[cut-1] = '\0';
		start = pos + 1;
		strArray1[i] = seq1;
		i++;
	}
	char * seq1 = (char *) space((cut + 1) * sizeof (char));	
	sscanf(start, "%s", seq1);
	start = pos;
	strArray1[i] = seq1;
	
	*numOfRNA = i+1;
    free(line);
    return;
}

/* remove the & from a string for two sequences */
PRIVATE char *tokenize_one(char *line)
{
    char *pos, *copy;
    int cut = -1;

    copy = (char *) space(strlen(line) + 1);
    (void) sscanf(line, "%s", copy);
    pos = strchr(copy, '&');
    if (pos)
    {
        cut = (int) (pos - copy) + 1;
        if (cut >= strlen(copy)) cut = -1;
        if (strchr(pos + 1, '&')) nrerror("more than one cut-point in input");
        for (; *pos; pos++) *pos = *(pos + 1); /* splice out the & */
    }
    if (cut > -1)
    {
        if (cut_point == -1) cut_point = cut;
        else if (cut_point != cut)
        {
            fprintf(stderr, "cut_point = %d cut = %d\n", cut_point, cut);
            nrerror("Sequence and Structure have different cut points.");
        }
    }
    free(line);
    return copy;
}

int comp_nums(const int *num1, const int *num2)
{
    if (*num1 < *num2) return -1;
    if (*num1 == *num2) return 0;
    if (*num1 > *num2) return 1;
    return 0;
}

/* get the values for option -u, write them in u_vals*/
/* max. length of the unstructured region has to be <= len longer seq.!!!*/

/* u_vals[u_vals[0]] contains the largest -u value <= len longer seq. */
int get_u_values(char unstrs[], int **u_vals, int l1)
{
    int min, max, tmp, uc, count, uunstr;
    char *token, *cp;

    if ((*u_vals) != NULL) free((*u_vals));
    (*u_vals) = (int*) space(102 * sizeof (int));
    if (unstrs[0] != '\0' && strchr(unstrs, '-'))
    {/*range contains symbol "-"*/
        const char delimiters[] = " -";

        if (strchr(unstrs, ','))
            nrerror("option -u : enter either a range using \"-\" or a comma seperated list\n");
        cp = strdup(unstrs);
        token = strtok(cp, delimiters);
        min = atoi(token);
        token = strtok(NULL, delimiters);
        max = atoi(token);
        free(cp);
        if (min > max)
        {
            tmp = min;
            min = max;
            max = tmp;
        } else if (min == max)
        {
            nrerror("option -u : you enterd a range where min = max, use min < max to define a range");
        }
        if (max - min > 100)
        {
            fprintf(stderr, "only the first 100 length value are used\n");
        }

        (*u_vals)[0] = (max - min + 1) <= 100 ? (max - min + 1) : 100;
        uc = 0;
        max = max < min + 99 ? max : min + 99;
        for (tmp = min; tmp <= max; tmp++)
        {
            if (tmp <= l1)
            {
                (*u_vals)[++uc] = tmp;
                /* printf("%d,",tmp); */
            } else
            {
                fprintf(stderr, "option -u: length %d is longer than length of longer sequence. Only values <= length of longer sequence are allowed.\n", tmp);
                break;
            }
        }
        (*u_vals)[0] = uc;
        if (uc < 1) return (0);
        return (1);
        /* comma seperated list of values, symbol "," */
    } else if (unstrs[0] != '\0' && strchr(unstrs, ','))
    {
        const char delimiters[] = " ,";
        if (strchr(unstrs, '-'))
            nrerror("option -u : enter either a range using \"-\" or a comma seperated list\n");

        cp = strdup(unstrs);
        token = strtok(cp, delimiters);
        uc = 1;
        (*u_vals)[1] = atoi(token);
        while ((token = strtok(NULL, delimiters)) && uc < 20)
            (*u_vals)[++uc] = atoi(token);
        if ((token = strtok(NULL, delimiters)))
        {
            fprintf(stderr, "the first 20 length value are used\n");
        }
        free(cp);
        (*u_vals)[0] = 0;
        uunstr = (uc) <= 100 ? uc : 100;
        qsort((*u_vals), (uunstr + 1), sizeof (int), (void *) comp_nums);
        for (count = 0; count < uunstr + 1; count++)
        {
            if ((*u_vals)[count] > l1)
            {
                fprintf(stderr, "option -u: length %d is longer than length of longer sequence. Only values <= length of longer sequence are allowed.\n", (*u_vals)[count]);
                break;
            }
        }
        (*u_vals)[0] = count - 1;
        if ((count - 1) < 1) return (0);
        return (1);
    } else if (unstrs[0] != '\0')
    {
        (*u_vals)[0] = 1;
        uunstr = atoi(unstrs);
        if (uunstr > l1) return (0);
        (*u_vals)[1] = uunstr;
        return (1);
    } else
    { /* default value */
        (*u_vals)[0] = 1;
        if (default_u > l1)
        {
            (*u_vals)[1] = l1;
            fprintf(stderr, "option -u = %d exceeds length of longer sequence, %d. -u is set length of longer sequence.\n", default_u, l1);
        }
        (*u_vals)[1] = default_u;
    }
    return 1;
}

/* divide the constraints string in intermolecular constrains (inter)
   and intramolecular constrains within both sequences */

/* len1 is the length of the LONGER input seq ! */
void seperate_bp(char **inter, int len1, char **intra_l, char **intra_s)
{
    int i, j, len;
    short *pt = NULL;
    char *temp_inter, *pt_inter;

    len = strlen((*inter));
    /* printf("inter\n%s\n",(*inter)); */
    i = len + 1;
    temp_inter = (char*) space(sizeof (char) *i);/* to make a pair_table convert <|> to (|) */
    pt_inter = (char*) space(sizeof (char) *i);
    /* if shorter seq is first seq in constrained string, write the
       longer one as the first one */
    temp_inter[strlen((*inter))] = '\0';
    pt_inter[strlen((*inter))] = '\0';
    if (cut_point < len1)
    {
        /* write the constrain for the longer seq first */
        for (j = 0, i = cut_point - 1; i < len; i++, j++)
        {
            switch ((*inter)[i])
            {
                case '(':
                    temp_inter[j] = ')';
                    pt_inter[j] = ')';
                    break;
                case ')':
                    temp_inter[j] = '(';
                    pt_inter[j] = '(';
                    break;
                default:
                    temp_inter[j] = (*inter)[i];
                    pt_inter[j] = '.';
            }
        }
        /* then add the constrain for the shorter seq */
        for (i = 0; i < cut_point - 1; i++, j++)
        {
            switch ((*inter)[i])
            {
                case '(':
                    temp_inter[j] = ')';
                    pt_inter[j] = ')';
                    break;
                case ')':
                    temp_inter[j] = '(';
                    pt_inter[j] = '(';
                    break;
                default:
                    temp_inter[j] = (*inter)[i];
                    pt_inter[j] = '.';
            }
        }
        cut_point = len1 + 1;
        strcpy((*inter), temp_inter);
    } else
    {
        for (i = 0; i < strlen((*inter)); i++)
        {
            switch ((*inter)[i])
            {
                case '(':
                    pt_inter[i] = '(';
                    break;
                case ')':
                    pt_inter[i] = ')';
                    break;
                default:
                    pt_inter[i] = '.';
            }
        }
    }

    pt = make_pair_table(pt_inter);

    /* intramolecular structure in longer (_l) and shorter (_s) seq */
    (*intra_l) = (char*) space(sizeof (char) *(len1 + 1));
    (*intra_s) = (char*) space(sizeof (char) *(strlen((*inter)) - len1 + 2));
    (*intra_l)[len1] = '\0';
    (*intra_s)[strlen((*inter)) - len1 + 1] = '\0';
    /* now seperate intermolecular from intramolecular bp */
    for (i = 1; i <= pt[0]; i++)
    {
        if (pt[i] == 0)
        {
            temp_inter[i - 1] = (*inter)[i - 1];
            if (i < cut_point)
            {
                (*intra_l)[i - 1] = (*inter)[i - 1];
                if ((*inter)[i - 1] == '|')
                    (*intra_l)[i - 1] = '.';
            } else
            {
                (*intra_s)[i - cut_point] = (*inter)[i - 1];
                if ((*inter)[i - 1] == '|')
                    (*intra_s)[i - cut_point] = '.';
            }
        } else
        {
            if (i < cut_point)
            {
                /* intermolekular bp */
                if (pt[i] >= cut_point)
                {
                    temp_inter[i - 1] = (*inter)[i - 1];
                    (*intra_l)[i - 1] = '.';
                    (*intra_s)[pt[i] - cut_point] = '.';
                } else
                { /* intramolekular bp */
                    (*intra_l)[i - 1] = (*inter)[i - 1];
                    temp_inter[i - 1] = '.';
                }
            } else
            { /* i>=cut_point */
                /* intermolekular bp */
                if (pt[i] < cut_point)
                {
                    temp_inter[i - 1] = (*inter)[i - 1];
                    /* (*intra_s)[i-1] = '.'; */
                } else
                { /* intramolekular bp */
                    (*intra_s)[i - cut_point] = (*inter)[i - 1];
                    temp_inter[i - 1] = '.';
                }
            }
        }
    }

    /* printf("%s -1\n%s -2\n%s -3\n%s -4\n",(*inter),temp_inter,(*intra_l),(*intra_s)); */
    strcpy((*inter), temp_inter);
    free(temp_inter);
    free(pt_inter);
    free(pt);
}

PRIVATE void print_interaction(interact *Int, char *s1, char *s2, pu_contrib *p_c, pu_contrib *p_c2, int w, int incr3, int incr5)
{
    char *i_long, *i_short;
    int i, len, l_l, l_s, len1, end5, end3, i_min, j_min, l1, add_a, add_b, nix_up;
    double p_c_S;
    double G_min, Gi_min, Gul, G_sum, Gus, diff;
    duplexT mfe;
    char *struc;

    G_min = Int->Gikjl;
    Gi_min = Int->Gikjl_wo;
    len1 = Int->length;
    len = strlen(s1) + strlen(s2);

    /* use duplexfold() to fold the interaction site */
    l_l = (Int->i - Int->k + 1);
    i_long = (char*) space(sizeof (char) *(l_l + 1));
    l_s = (Int->l - Int->j + 1);
    i_short = (char*) space(sizeof (char) *(l_s + 1));

    strncpy(i_long, &s1[Int->k - 1], l_l);
    i_long[l_l] = '\0';
    strncpy(i_short, &s2[Int->j - 1], l_s);
    i_short[l_s] = '\0';

    mfe = duplexfold(i_long, i_short);

    i_min = mfe.i;
    j_min = mfe.j;
    l1 = strchr(mfe.structure, '&') - mfe.structure;

    /* printf("%s %3d,%-3d : %3d,%-3d (%5.2f)\n", mfe.structure, i_min+1-l1,
       i_min, j_min, j_min+strlen(mfe.structure)-l1-2, mfe.energy ); */

    /* structure by duplexfold is shorter than structure by RNAup:*/

    add_a = add_b = 0; /* length difference in longer / shorter sequence*/
    nix_up = 0;
    if (((i_min + 1 - l1) - i_min) != (Int->k - Int->i))
    {
        add_a = Int->i - Int->k + 2;
    }
    if (((j_min + strlen(mfe.structure) - l1 - 2) - j_min) != (Int->l - Int->j))
    {
        add_b = Int->l - Int->j + 2;
    }
    /* printf("add_a %d   add_b %d\n",add_a,add_b); */
    if (add_a || add_b)
    {
        nix_up = 1;
        if (add_a && add_b == 0) add_b = Int->l - Int->j + 2;
        if (add_a == 0 && add_b) add_a = Int->i - Int->k + 2;
        struc = (char*) space(sizeof (char) *(add_a + add_b + 3));
        for (i = 0; i < (add_a + add_b - 1); i++)
        {
            if (i != l_l) struc[i] = '.';
            if (i == l_l) struc[i] = '&';
        }
        struc[i] = '\0';
    } else
    {
        l1 = strlen(mfe.structure);
        struc = (char*) space(sizeof (char) *(l1 + 1));
        strcpy(struc, mfe.structure);
    }

    end5 = MAX(1, Int->k - incr5);
    end3 = MIN(MIN(l_l - 1 + incr3, w + incr3 + incr5), len1);
    p_c_S = p_c->H[end5][end3] + p_c->I[end5][end3] + p_c->M[end5][end3] + p_c->E[end5][end3];
    Gul = -RT * log(p_c_S);

    if (p_c2 == NULL)
    {
        G_sum = Gi_min + Gul;

        /* printf("dG = dGint + dGu_l\n"); */
        // printf("%s %3d,%-3d : %3d,%-3d (%.2f = %.2f + %.2f)\n",
        //         struc, Int->k, Int->i, Int->j, Int->l, G_min, Gi_min, Gul);
        // printf("%s&%s\n", i_long, i_short);
    } else
    {
        p_c_S = p_c2->H[Int->j][(Int->l)-(Int->j)] +
                p_c2->I[Int->j][(Int->l)-(Int->j)] +
                p_c2->M[Int->j][(Int->l)-(Int->j)] +
                p_c2->E[Int->j][(Int->l)-(Int->j)];
        Gus = -RT * log(p_c_S);
        G_sum = Gi_min + Gul + Gus;
        /* printf("dG = dGint + dGu_l + dGu_s\n"); */
        // printf("%s %3d,%-3d : %3d,%-3d (%.2f = %.2f + %.2f + %.2f)\n",
        //         struc, Int->k, Int->i, Int->j, Int->l, G_min, Gi_min, Gul, Gus);
        // printf("%s&%s\n", i_long, i_short);
    }
    if (!EQUAL(G_min, G_sum))
    {
        printf("ERROR\n");
        diff = fabs((G_min)-(G_sum));
        printf("diff %.18f\n", diff);
    }
    if (nix_up) fprintf(stderr, "RNAduplex structure doesn't match any structure of RNAup structure ensemble\n");
    free(i_long);
    free(i_short);
    free(mfe.structure);
    free(struc);
}

/* print coordinates and free energy for the region of highest accessibility */
PRIVATE void print_unstru(pu_contrib *p_c, int w)
{
    int i, j, len, min_i, min_j;
    double dG_u, min_gu;

    if (p_c != NULL)
    {
        min_gu = 1000.0;
        len = p_c->length;

        for (i = 1; i <= len; i++)
        {
            for (j = i; j < MIN((i + w), len + 1); j++)
            {
                double blubb;
                if ((j - i + 1) == w && i + w - 1 <= len)
                {
                    blubb = p_c->H[i][j - i] + p_c->I[i][j - i] + p_c->M[i][j - i] + p_c->E[i][j - i];
                    dG_u = -RT * log(blubb);
                    if (dG_u < min_gu)
                    {
                        min_gu = dG_u;
                        min_i = i;
                        min_j = i + w - 1;
                    }
                }
            }
        }
        printf("%d,%d \t (%.3f) \t for u=%d\n", min_i, min_j, min_gu, w);
    } else
    {
        nrerror("error with prob unpaired");
    }
}


void reverseRNAupOutput(int location, int n1, int n2, double ***** rnaupCollections)	//rnaWins
{

	int i,j,w;
//printf("** A \n");

	/*double *** matrix = malloc(sizeof(double **) * MAX_RNA_SIZE);
	for(i=1;i<=n1;i++)
	{
		matrix[i] = malloc(sizeof(double *) * MAX_RNA_SIZE);
		for(j=1;j<=n2;j++)
		{
			matrix[i][j] = malloc(sizeof(double) * 26);
		}
	}
*/	
//printf("** B \n");
	for(i=1;i<=n1;i++)
	{
		for(j=n2;j>=1;j--)
		{
			
			for(w=1;w<=25;w++)
			{
				int w2;
				for(w2=1;w2<=25;w2++)
				{
					//printf("i=%d, j=%d, w=%d \n", i,j,w);
					int pseudo_L = n2 - j + 1;
					int pseudo_R = pseudo_L + w2;
					
					if(pseudo_R > n2)
						continue;
					//else printf("X \n");
					matrix[i][pseudo_R][w][w2] = rnaupCollections[location][i][j][w][w2];
				}
			}
		}
	}
//printf("** C \n");	
	for(i=1;i<=n1;i++)
	{
		for(j=1;j<=n2;j++)
		{
			for(w=1;w<=25;w++)
			{
				int w2;
				for(w2=1;w2<=25;w2++)
				{
					rnaupCollections[location][i][j][w][w2] = matrix[i][j][w][w2];
				}
			}
		}
	}
//printf("** D \n");
/*
	for(i=1;i<=n1;i++)
		for(j=1;j<=n2;j++)
			if(matrix[i][j] != NULL)
				free(matrix[i][j]);

	for(i=1;i<=n1;i++)
		if(matrix[i] != NULL)
			free(matrix[i]);

	if(matrix != NULL)
		free(matrix);
*/
}



