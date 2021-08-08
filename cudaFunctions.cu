#include <cuda_runtime.h>
#include <helper_cuda.h>
#include "myProto.h"

__device__ int myStrlen(char arr[]){ 
	int counter=0;
	for(int i=0;i<7;i++){// max size in Semiconservative group is 5
		if(arr[i]=='\0')
			return counter;
		else
			counter++;
	}
	return counter;
}

//Computing values and overrides vec A
__global__ void findMuTantCuda(char *d_A, char *seq1,int offset,int isMax,int seq2Size,char *d_Seq2Temp,int *d_weigh)
{
    int index = blockDim.x * blockIdx.x + threadIdx.x;
    if(index<seq2Size){// need only strlen(size2) trheads to do the job
    	if(isMax==1){ //its maximum
    		int flag1 = 0, flag2 = 0, i, j;
			char ConservativeGroup[][5] = { { 'S','T','A','\0' } ,{ 'N','E','Q','K','\0' } ,{ 'N','D','E','Q','\0'},
			{ 'N','H','Q','K','\0' },{ 'Q','H','R','K','\0' },{ 'M','I','L','V','\0' },{ 'M','I','L','F','\0' },{ 'H','Y','\0' },{ 'F','Y','W','\0'} };
		
			for (i = 0; i < 9; i++) { //checking if the letters are in the same conservative group
				int size = myStrlen(ConservativeGroup[i]);// checking the size of the next group of chars
				for (j = 0; j < size; j++) {
			
					if (seq1[index+offset] == ConservativeGroup[i][j])//first char found in a group of close attributes
						flag1 = 1;
					if (d_A[index] == ConservativeGroup[i][j])//second char found in a group of close attributes
						flag2 = 1;
					if (flag1 == 1 && flag2 == 1)// the chars are in the same conservative group
					{
						d_Seq2Temp[index]=d_A[index];
						return;
					}
				}
				flag1 = 0;
				flag2 = 0;
				size = 0;
			}
			d_Seq2Temp[index]=seq1[index+offset]; //if not in the same conservative group return the letter from seq1
    	
    	 }
    	 else{//its minimum
    	 	if(d_weigh[2]>=d_weigh[3]){ //try to get the max weight to get minimum result
    	 	    char semiConservativeGroup[11][7] = { { 'C','S','A','\0' } ,{ 'A','T','V','\0' } ,{ 'S','A','G','\0'},
				{ 'S','T','N','K','\0'},{ 'S','T','P','A','\0'},{ 'S','G','N','D','\0'},{ 'S','N','D','E','Q','K','\0' },{ 'N','D','E','Q','H','K','\0'}
				,{'N','E','Q','H','R','K','\0'},{'F','V','L','I','M','\0' },{'H','F','Y','\0' } };
				for (int i = 0; i < 9; i++) { //checking if the letters are in the same Semiconservative group
					int size = myStrlen(semiConservativeGroup[i]);// checking the size of the next group of chars
					for (int j = 0; j < size; j++) {
						if (seq1[index+offset] == semiConservativeGroup[i][j])//first char found in a Semiconservative group
						{
							if(j==0)
								d_Seq2Temp[index]= semiConservativeGroup[i][1];
							else
								d_Seq2Temp[index]= semiConservativeGroup[i][0];
							return;
						}

					}
					size = 0;
				}
				
    	 	 }

			if(seq1[index+offset]!='B') // if the char cant be fount in any Semiconservative group then return letter that cant be found in any group
				d_Seq2Temp[index]= 'B';
			else
				d_Seq2Temp[index]= 'J';
		
    	 }
    	 
    }

}

int computeOnGPU(char *seq1,char *seq2,int numOfElementSeq1,int numOfElementsSeq2,int offset,int isMax,int *weight)
{

    // Error code to check return values for CUDA calls
    cudaError_t err = cudaSuccess;
    size_t sizeSeq1 =numOfElementSeq1 * sizeof(char);
    size_t sizeSeq2 = numOfElementsSeq2 * sizeof(char);
    size_t sizeWeight=4*sizeof(int);
  


    char *d_Seq2;
    char *d_Seq1;
    char *d_Seq2Temp;
    int *d_weight;
    err = cudaMalloc((void **)&d_Seq2, sizeSeq2);     // Allocate memory on GPU to copy seq1 from the host
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    err = cudaMalloc((void **)&d_Seq1, sizeSeq1);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
        err = cudaMalloc((void **)&d_Seq2Temp, sizeSeq2);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
        err = cudaMalloc((void **)&d_weight, sizeWeight);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to allocate device memory - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy data from host to the GPU memory
    err = cudaMemcpy(d_Seq2, seq2, sizeSeq2, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
        // Copy data from host to the GPU memory
    err = cudaMemcpy(d_Seq1, seq1, sizeSeq1, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
        err = cudaMemcpy(d_weight, weight, sizeWeight, cudaMemcpyHostToDevice);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy data from host to device - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
    
    // Launch the Kernel
    int threadsPerBlock = 256;
    int blocksPerGrid = (numOfElementsSeq2 + threadsPerBlock - 1) / threadsPerBlock;
    findMuTantCuda<<<blocksPerGrid, threadsPerBlock>>>(d_Seq2, d_Seq1,offset,isMax,numOfElementsSeq2,d_Seq2Temp,d_weight);
    
    
    err = cudaGetLastError();
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to launch vectorAdd kernel -  %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Copy the  result from GPU to the host memory.
    err = cudaMemcpy(seq2, d_Seq2Temp, sizeSeq2, cudaMemcpyDeviceToHost);
    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy result array from device to host -%s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    // Free allocated memory on GPU
    if (cudaFree(d_Seq2) != cudaSuccess)
    {
       fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
        // Free allocated memory on GPU
    if (cudaFree(d_Seq2Temp) != cudaSuccess)
    {
       fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
        // Free allocated memory on GPU
    if (cudaFree(d_Seq1) != cudaSuccess)
    {
       fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
        if (cudaFree(d_weight) != cudaSuccess)
    {
       fprintf(stderr, "Failed to free device data - %s\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    return 0;
}
