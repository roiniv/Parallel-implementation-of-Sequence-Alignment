
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include "myProto.h"
#include <mpi.h>
#define SEQ1 10000
#define SEQ2 5000

MPI_Datatype bestMutant_type;

/* create our bestMutant struct type*/
void createBestMutantType(MPI_Datatype *bestMutant_type)
{
	int block_lengths[3] = {1, SEQ2, 1};
	MPI_Aint disp[3];
	MPI_Datatype types[3] = {MPI_INT,MPI_CHAR,MPI_INT};

	disp[0] = offsetof(BestMutant,offset);
	disp[1] = offsetof(BestMutant,mutant);
	disp[2] = offsetof(BestMutant,score);


	MPI_Type_create_struct(3,block_lengths,disp,types,bestMutant_type);
	MPI_Type_commit(bestMutant_type);
}


int main(int argc, char *argv[])
{
    int rank, i, size;
    int *weight=(int *)malloc(sizeof(int)*4);
    int *isMax=(int*)malloc(sizeof(int)); //==1 if max and 0 if min
    char *seq1=(char*)malloc(sizeof(char)*SEQ1);
    char *seq2=(char*)malloc(sizeof(char)*SEQ2);
    char *seq2Temp;
    int startingOffset=0, endingOffset=0,numOfOffset;;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) //can run with only two processes!
    {
       printf("Run this project with two processes only\n");
       MPI_Abort(MPI_COMM_WORLD, __LINE__);
    }

    if (seq2 == NULL||seq1 == NULL||isMax == NULL||weight == NULL) //check if all malloc Successed
    	MPI_Abort(MPI_COMM_WORLD, __LINE__);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    char filename[25] = "input.txt";

    if (rank == 0){  // Reading from file and divide the tasks between both processes

        readFromFile(filename,weight ,seq1,seq2,isMax);
        numOfOffset=strlen(seq1)-strlen(seq2)+1;
        int start=numOfOffset/2; //set the start offset for the second process
        MPI_Send(&start, 1, MPI_INT, 1, 0, MPI_COMM_WORLD); //Send half of the offsets to the  second process
        endingOffset=start; //set the ending offset for the first process ass the start offset of the second process

    }
    MPI_Bcast(seq1, SEQ1, MPI_CHAR, 0, MPI_COMM_WORLD); //Broadcasting seq1
    MPI_Bcast(seq2, SEQ2, MPI_CHAR, 0, MPI_COMM_WORLD); //Broadcasting seq2
    MPI_Bcast(weight, 4, MPI_INT, 0, MPI_COMM_WORLD); //Broadcasting weight array
    MPI_Bcast(isMax, 1, MPI_INT, 0, MPI_COMM_WORLD); //Broadcasting isMax
    MPI_Bcast(&numOfOffset, 1, MPI_INT, 0, MPI_COMM_WORLD); //Broadcasting the number of offsets

    if(rank!=0){
    	 MPI_Recv(&startingOffset, 1 , MPI_INT, 0, 0, MPI_COMM_WORLD, &status); //recv half of the offsets to work on
    	 endingOffset=numOfOffset;
    }

    seq2Temp=(char*)malloc(sizeof(char)*strlen(seq2)); //create temp from seq2 so we dont change seq2 directly

    BestMutant *bestMutant; //this struct will hold the best mutant result
    bestMutant=(BestMutant*)malloc(sizeof(BestMutant));
    bestMutant->offset=-1;

    if (seq2Temp == NULL||bestMutant==NULL)
    	MPI_Abort(MPI_COMM_WORLD, __LINE__);

    createBestMutantType(&bestMutant_type);// create this struct type

    for(int offset=startingOffset;offset<endingOffset;offset++){ //this loop run, for each process ,half of the number of offsets
    	strcpy(seq2Temp,seq2);
    	#pragma omp parallel for //for each offset, get half of the mutant using openMp
    	   for (i = 0; i <strlen(seq2)/2; i++)
    		   findMuTantOmp(seq1,seq2Temp,offset,i,*isMax,weight);

    	   /*for each offset, get second half of the mutant using Cuda-
    	    *  if the size of seq2 is even send strlen(seq1)-strlen(seq2)/2 elemnts else send one more*/
    	   if(strlen(seq2)%2==0) {
    		if (computeOnGPU(seq1+strlen(seq2)/2,seq2Temp + strlen(seq2)/2,strlen(seq1)-strlen(seq2)/2,strlen(seq2)/2 ,offset, *isMax,weight) != 0)
    			MPI_Abort(MPI_COMM_WORLD, __LINE__);
    	}
    	else{
    		if (computeOnGPU(seq1+strlen(seq2)/2,seq2Temp + strlen(seq2)/2,strlen(seq1)-strlen(seq2)/2+1, strlen(seq2)/2+1,offset, *isMax,weight) != 0)
    			MPI_Abort(MPI_COMM_WORLD, __LINE__);
    	}

    	int score=getScore(seq1,seq2Temp,offset,weight); // calculate the score of each mutant
    	if((bestMutant->offset==-1)||(*isMax==1&&score>bestMutant->score)||(*isMax==0&&score<bestMutant->score)){ //save the best result in the bestMustant struct
        	strcpy(bestMutant->mutant,seq2Temp);
        	bestMutant->offset=offset;
        	bestMutant->score=score;
    	}
    }

    // Collect the result on one of processes and write it to the output file
    if (rank == 0){
    	BestMutant bestMutantFrom1;
    	MPI_Recv(&bestMutantFrom1, 1, bestMutant_type, 1, 0, MPI_COMM_WORLD, &status); //recv best result from the second process
    	if((bestMutantFrom1.score>bestMutant->score&&*isMax==1)||(bestMutantFrom1.score<bestMutant->score&&*isMax==0)){ //save the best result
    		bestMutant->score=bestMutantFrom1.score;
    		bestMutant->offset=bestMutantFrom1.offset;
    		strcpy(bestMutant->mutant,bestMutantFrom1.mutant);

    	}
    	strcpy(filename ,"output.txt");
    	writeToFile(bestMutant, filename); //write the result to the output file

    }

    else
        MPI_Send(bestMutant, 1, bestMutant_type, 0, 0, MPI_COMM_WORLD); //send the best result to the first process



    MPI_Finalize();
    return 0;
}
