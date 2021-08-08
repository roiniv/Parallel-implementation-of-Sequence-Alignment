#pragma once


typedef struct BestMutant{
	int offset;
	char mutant[5000];
	int score;

}BestMutant;


int computeOnGPU(char *seq1,char *seq2,int numOfElementSeq1,int numOfElementsSeq2,int offset,int isMax,int *weight);
void readFromFile(char *filename, int *weight,char *seq1,char *seq2,int *isMax);
int isSemiConservativeGroup(char firstC, char secondC);
int getScore(char *seq1,char *seq2,int offset,int *weight);
int isConservativeGroup(char firstC, char secondC);
void findMuTantOmp(char *seq1,char *seq2, int offset,int index, int isMax,int *weight);
char changeMaxLetter(char a,char b);
char changeMinLetter(char a,int *weight);
char getSemiConservativeGroupletter(char c);
void writeToFile(BestMutant *bestMutant, char * filename);

