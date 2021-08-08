#include "myProto.h"
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#define SEQ1 10000
#define SEQ2 5000




//return 1 if firstC and secondC is in the same Semi conservative group; other return 0
int isSemiConservativeGroup(char firstC, char secondC)//checks if two chars have a similar attributes
{
	int flag1 = 0, flag2 = 0, i, j;
	char close[11][7] = { { 'C','S','A','\0' } ,{ 'A','T','V','\0' } ,{ 'S','A','G','\0'},
	{ 'S','T','N','K','\0'},{ 'S','T','P','A','\0'},{ 'S','G','N','D','\0'},{ 'S','N','D','E','Q','K','\0' },{ 'N','D','E','Q','H','K','\0'}
	,{'N','E','Q','H','R','K','\0'},{'F','V','L','I','M','\0' },{'H','F','Y','\0' } };
	for (i = 0; i < 11; i++) {
		int size = strlen(close[i]);// checking the size of the next group of chars
		for (j = 0; j < size; j++) {
			if (firstC == close[i][j]) //first char found in a group of similar attributes
				flag1 = 1;
			if (secondC == close[i][j])//second char found in a group of similar attributes
				flag2 = 1;
			if (flag1 == 1 && flag2 == 1) // the chars are similiar
				return 1;
		}

		flag1 = 0;
		flag2 = 0;
	}
	return 0; //the chars not similar
}

//return 1 if firstC and secondC is in the same conservative group; other return 0
int isConservativeGroup(char firstC, char secondC)
{
	int flag1 = 0, flag2 = 0, i, j;
	char close[][5] = { { 'S','T','A','\0' } ,{ 'N','E','Q','K','\0' } ,{ 'N','D','E','Q','\0'},
	{ 'N','H','Q','K','\0' },{ 'Q','H','R','K','\0' },{ 'M','I','L','V','\0' },{ 'M','I','L','F','\0' },{ 'H','Y','\0' },{ 'F','Y','W','\0'} };
	for (i = 0; i < 9; i++) {

		int size = strlen(close[i]);// checking the size of the next group of chars
		for (j = 0; j < size; j++) {
			if (firstC == close[i][j])//first char found in a group of close attributes
				flag1 = 1;
			if (secondC == close[i][j])//second char found in a group of close attributes
				flag2 = 1;
			if (flag1 == 1 && flag2 == 1)// the chars have close attributes
				return 1;
		}

		flag1 = 0;
		flag2 = 0;
		size = 0;
	}
	return 0;// the chars not  have close attributes

}


char getSemiConservativeGroupletter(char c)//search for letter that is in the same group as c ,return null if their is no such letter
{
	int  i, j;
	char close[11][7] = { { 'C','S','A','\0' } ,{ 'A','T','V','\0' } ,{ 'S','A','G','\0'},
	{ 'S','T','N','K','\0'},{ 'S','T','P','A','\0'},{ 'S','G','N','D','\0'},{ 'S','N','D','E','Q','K','\0' },{ 'N','D','E','Q','H','K','\0'}
	,{'N','E','Q','H','R','K','\0'},{'F','V','L','I','M','\0' },{'H','F','Y','\0' } };
	for (i = 0; i < 9; i++) {
		int size = strlen(close[i]);// checking the size of the next group of chars
		for (j = 0; j < size; j++) {
			if (c == close[i][j])//first char found in a group of close attributes
			{
				if(j==0)
					return close[i][1];
				else
					return close[i][0];
			}

		}
		size = 0;
	}
	return NULL;

}

/* this func changing only one letter in seq2 . it find the best letter according to the rules*/
void findMuTantOmp(char *seq1,char *seq2, int offset,int index, int isMax,int *weight)
{
    if(isMax==1){
    	seq2[index]=changeMaxLetter(seq1[offset+index],seq2[index]);
    }
    else{
    	seq2[index]=changeMinLetter(seq1[offset+index],weight);
    }

}


/* by two giving chars returns the minimum letter according to the rules,
 /*if w3>w4 then we need to find letter that in the same semi conservative group as char a(seq1) else we need to return letter that cant be found in any group */
char changeMinLetter(char a,int *weight){
	if(weight[2]>=weight[3]){

		char c=getSemiConservativeGroupletter(a);
		if(c!=NULL)
			return c;
	}

	if(a!='B')
		return 'B';
	else
		return 'J';
}

/* by two giving chars return the maximum letter according to the rules,
 * first check is the letter in the same conservative group, if they are so not change will be accepted,other change for the same latter as the first one*/
char changeMaxLetter(char a,char b){
	if(isConservativeGroup(a,b)==1)
		return b;
	else
		return a;
}

// by two giving  sequencers(seq1,seq2),the offset and weight array this func return the alignment score(using openMP)
int getScore(char *seq1,char *seq2,int offset,int *weight){
	int w1Counter=0,w2Counter=0,w3Counter=0,w4Counter=0;

	//#pragma omp parallel for //using openMp
	for (int i = 0; i < strlen(seq2); i++){
		if(seq1[i+offset]==seq2[i])
			w1Counter++;
		else if(isConservativeGroup(seq1[i+offset],seq2[i])==1)
			w2Counter++;
		else if(isSemiConservativeGroup(seq1[i+offset],seq2[i])==1)
			w3Counter++;
		else
			w4Counter++;

	}
	return (w1Counter*weight[0])-(w2Counter*weight[1])-(w3Counter*weight[2])-(w4Counter*weight[3]);
}


/*read from file the weights, two sequences and maximum/minimum */
void readFromFile(char *filename, int *weight,char *seq1,char *seq2,int *isMax) //Read from file function
{
    double *temp;
    char seq1Temp[SEQ1];
    char seq2Temp[SEQ2];
    char maxOrMin[10];
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "failed to open %s", filename);
        exit(0);
    }
    fscanf(fp, "%d %d %d %d", &weight[0], &weight[1], &weight[2], &weight[3]);
    fscanf(fp, "%s", seq1Temp);
    fscanf(fp, "%s", seq2Temp);
    fscanf(fp, "%s", maxOrMin);
    seq1=(char*)realloc(seq1,sizeof(char)*strlen(seq1Temp));
    seq2=(char*)realloc(seq2,sizeof(char)*strlen(seq2Temp));
    strcpy(seq1,seq1Temp);
    strcpy(seq2,seq2Temp);
    if(strcmp(maxOrMin,"maximum")==0)
    	*isMax=1;
    else
    	*isMax=0;
    fclose(fp);

}

/*write the result to the file, first lain contains best mutant,secnod line contains is offset and score*/
void writeToFile(BestMutant *bestMutant, char * filename){

	FILE* f = fopen(filename,"w");
	if(!f){
		printf("cant open file");
		exit(-1);
	}
	fprintf(f, "Best mutant is : %s\n", bestMutant->mutant);
	fprintf(f,"The offset is :%d, and the score is : %d \n", bestMutant->offset,bestMutant->score);
	fclose(f);

}
