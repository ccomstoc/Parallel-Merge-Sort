//Maddie Neely and Connor Comstock
//Project 3

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <bits/stdc++.h>
#include <cmath>
#include "mpi.h" // message passing interface

using namespace std;

// New compile and run commands for MPI!
// mpicxx -o blah file.cpp
// mpirun -q -np 32 blah

//helpful functions for filling array and printing arrays neatly for debugging
void fill(int* a, int size);
void print(string str, int id, int* a, int n);


//Here are the functions we will need to write...
void mergesort (int * a, int first, int last);
void smerge(int * a, int * b, int lasta, int lastb, int * output = NULL); //wrote in last project
int rankF(int * a, int first, int last, int valToFind);
void pmerge(int * a, int * b, int lasta, int lastb, int * output = NULL);


int my_rank;			// my CPU number for this process
int p;					// number of CPUs that we have

int main (int argc, char * argv[]) {


	int source;				// rank of the sender
	int dest;				// rank of destination
	int tag = 0;			// message number
	char message[100];		// message itself
	MPI_Status status;		// return status for receive

	// Start MPI
	MPI_Init(&argc, &argv);

	// Find out my rank!
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	// Find out the number of processes!
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	// THE REAL PROGRAM IS HERE

	int seed = time(0);
	srand(seed);

	//int bookExample [] = {4,6,8,9,16,17,18,19,20,21,23,25,27,29,31,32,1,2,3,5,7,10,11,12,13,14,15,22,24,26,28,30};
  //input = {1,3,4,6,7,9,10,11,13,15,16,17,20,21,22,28,2,5,8,12,14,18,19,23,24,25,26,27,29,30,31,32};

	int n2 = 0;//holds user inputed array size

	//first populate array with some random data.... (integers) on process 0...
	if (my_rank == 0) // have process 0 create initial array and broadcast it to every other process
	{
		cout << "Please enter the size of the array as a power of 2: " << endl;
		cin >> n2;

		while (ceil(log2(n2)) != floor(log2(n2)))//ensures input is a power of 2
		{
			cout << "Please enter the size of the array as a power of 2: " << endl;
			cin >> n2;
			cout << " n2 " << n2;

		}


	}
	MPI_Bcast(&n2, 1, MPI_INT, 0, MPI_COMM_WORLD);


	int * sad = new int[n2];
	if (my_rank == 0) // have process 0 create initial array and broadcast it to every other process
	{
		fill(sad,n2);
	}
	MPI_Bcast(sad, n2, MPI_INT, 0, MPI_COMM_WORLD);

	if(my_rank == 0)
		print("Before Sort: ", my_rank, sad, n2);

	mergesort(sad,0,n2-1);

	if(my_rank == 0)
		print("AFTER MergeSort", my_rank, sad, n2);

  delete [] sad; //yay no mem leks

	// Shut down MPI
	MPI_Finalize();

	return 0;
}

void fill(int* a, int size) {//fill an array with randoms integers up to 500
	int i;
	for (i = 0; i < size; i++) {
		a[i] = rand() % 500; //;
	}
}

void print(string str, int id, int* a, int n) {
	cout << str << ", Process " << id << " has: " << endl;
	for(int i=0; i < n; i++)
	{
		cout << a[i] << " ";
	}
	cout << endl;
	cout << endl;
}



void mergesort (int * a, int first, int last){

	if(last - first == 0) //condition to escape recursion
		return;

	int middle = ((last-first)/2)+first; //find position of middle
	int * output = new int[last+1]; //create output array

	mergesort(a, first, middle); //call mergesort on first half
	mergesort(a, middle+1, last); //call on second half

	pmerge(&a[first], &a[middle+1], middle-first, last-middle-1, output); //do the merge


	for(int i = first; i <= last; i++)
	{
		a[i] = output[i-first]; //copy output back to a
	}

	delete [] output; //delete output because we are not heathens
}

//Serial MergeSort
void smerge(int * a, int * b, int lasta, int lastb, int * output){

	int i,j,k; //create and set iterators to 0
	i = j = k = 0;

	while ((i <= lasta) && (j <= lastb)) //not at end of either "array"
	{
		if(a[i] < b[j]) //compare and add the lower value first
		{
			output[k++] = a[i++];
		}
		else
		{
			output[k++] = b[j++];
		}
	}
	//check which array didn't reach the end and then add the remaining values
	while(i <= lasta)
	{
		output[k++] = a[i++];
	}
	while(j <= lastb)
	{
		output[k++] = b[j++];
	}
}

int rankF(int * a, int first, int last, int valToFind) //a recursive function that is basically binary search
{
	/*Rank function determines where "valToFind" would be potioned in array "a". first and last dictate
	the bounds to search*/
	//Basically binary search
	 if (last >= first) {
        int middle = ((last-first)/2)+first; //find position of middle

         if (last-first == 0)
         {
         	if(valToFind > a[last])
         	{
         		return last+1;
         	}
         	return last;
         }
        // If element is smaller than mid then it's in left subarrray
        if (a[middle] > valToFind)
            return rankF(a, first, middle, valToFind);

        // Otherwise it must be in right subarray
        return rankF(a, middle + 1, last, valToFind);
    }

    // Just in case element is not present...
    return -1;
}

void pmerge(int * a, int * b, int lasta, int lastb, int * output)
{
	//When there is only one element in each
	if(lasta == 0 && lastb==0)
	{
		if(a[lasta]>b[lastb])
		{
			output[0] = b[lastb];
			output[1] = a[lasta];
		}
		else
		{
			output[1] = b[lastb];
			output[0] = a[lasta];
		}
		return;
	}

	//----------calculate rank in parallel----------
	double stepSizeA;
	double stepSizeB;
	stepSizeA = log2(double(lasta+1));
	stepSizeB = log2(double(lastb+1));

	int * selectA = new int[int(stepSizeA)];
	int * selectB = new int[int(stepSizeB)];

	int * SRANKA = new int[int(stepSizeA)];
	int * SRANKB = new int[int(stepSizeB)];

	//step 1 is to caluclate SRANKA and SRANKB //we assume input is base 2

	for(int i = 0; i<int(stepSizeA); i++)
	{
		selectA[i] = a[i*int(stepSizeA)];
		selectB[i] = b[i*int(stepSizeB)];
	}
	//WIll store local info that we will smash together and share
	int * localSRANKA = new int[int(stepSizeA)] ();
	int * localSRANKB = new int[int(stepSizeA)] ();

	int batchNum = 0;//used to keep track of parallel iterations

	//uses *math* to determine which rank this processer needs to calculate
	while((my_rank+batchNum) < int(stepSizeA))
	{

		localSRANKA[my_rank+batchNum] = rankF(b,0,lastb,selectA[my_rank+batchNum]);
		localSRANKB[my_rank+batchNum] = rankF(a,0,lasta,selectB[my_rank+batchNum]);

		batchNum += p;
	}
	//Smashes the srank arrays into one that is shared with all processors
	MPI_Allreduce(localSRANKA,SRANKA,stepSizeA,MPI_INT, MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(localSRANKB,SRANKB,stepSizeB,MPI_INT, MPI_MAX,MPI_COMM_WORLD);


	//----------Set up win array and final indicies----------
	int * WIN = new int[lasta+lastb+2] (); //final mergered array, +2 bc indexed at 0
	int * Windex = new int[int(stepSizeA)*2] ();//Index of final local array
	int * FinalWindex = new int[int(stepSizeA)*2] ();//Index of final global array
	int * localWin = new int[lasta+lastb+2] ();//Local processers part of win array

	for(int x = 0; x<int(stepSizeA); x++)
	{
		int aindex = x*int(stepSizeA)+SRANKA[x];//Calculate postion it will be in win array
		int bindex = x*int(stepSizeA)+SRANKB[x];//Does it simotaniously for a and b
		WIN[aindex] = a[x*int(stepSizeA)];
		WIN[bindex] = b[x*int(stepSizeB)];

		if(my_rank == 0)
		{
			localWin[aindex] = a[x*int(stepSizeA)];
			localWin[bindex] = b[x*int(stepSizeB)];
		}

		Windex[x] = aindex;
		Windex[x+int(stepSizeA)] = bindex;

	}
	//Could replace with pmerge when implimentation is finished, though time added in minimal
	smerge(&Windex[0],&Windex[int(stepSizeA)],int(stepSizeA)-1,int(stepSizeA)-1,FinalWindex);

	//----------SHAPES----------

	//Start by collecting all of the boundries in arrays
	int * shapeCordA = new int[2*int(stepSizeA)*2];
	int * shapeCordB = new int[2*int(stepSizeA)*2];

	//Add SelectPos and SelectPos+1 To Both Cord Arrays, Select is never apart of shapes, the first will act as a floor, and then + 1 for the next ceiling
	for(int i = 0; i < (int)stepSizeA*2; i+=2){
		int stepSizeElement = ((i/2)*stepSizeA);//((i/2)*stepSizeA), goes thru each selectPos element, in this example 0,4,8,12

		if(stepSizeElement != 0){//Skip zero because it acts as a floor for a shape that never exsists
			shapeCordA[i-1] = stepSizeElement;
			shapeCordA[i] = stepSizeElement+1;

			shapeCordB[i-1] = stepSizeElement;
			shapeCordB[i] = stepSizeElement+1;
		}else{
			shapeCordA[i] = 1;//we always skip zero but we want zero +1, it the ceiling for the first shape
			shapeCordB[i] = 1;
		}

		if(i == ((int)stepSizeA*2)-2 ){//When on the last element, add a floor for the last shape
			shapeCordA[i+1] = (lasta+1);
			shapeCordB[i+1] = (lasta+1);
		}
	}

	//Add according Rank Positions twice To Both Cord Arrays(They act as floor and ceiling, so one for each)
	for(int i = (int)stepSizeA*2; i < ((int)stepSizeA*2)*2; i+=2){
		//cout << "index: " << i << endl;
		int sRankElement = (i-(int)stepSizeA*2)/2;//for iterating thru srank, in this case gives, 0-3
		shapeCordA[i] = SRANKB[sRankElement];//
		shapeCordA[i+1] = SRANKB[sRankElement];

		shapeCordB[i] = SRANKA[sRankElement];//
		shapeCordB[i+1] = SRANKA[sRankElement];
	}

	//sort CordArrays(We can use merge becasue we built the arrays with both halfs sorted)
	int * mergeArray = new int[2*2*int(stepSizeA)];
	int numShapeCoord = 2*2*int(stepSizeA);

	//Smerge B
	//Could replace with pmerge when implimentation is finished
	smerge(&shapeCordB[0],&shapeCordB[numShapeCoord/2],(numShapeCoord/2)-1,(numShapeCoord/2)-1,mergeArray);

	for(int i = 0; i < numShapeCoord; i++)
		shapeCordB[i] = mergeArray[i];

	//Smerge A
	smerge(&shapeCordA[0],&shapeCordA[numShapeCoord/2],(numShapeCoord/2)-1,(numShapeCoord/2)-1,mergeArray);
		for(int i = 0; i < numShapeCoord; i++)
			shapeCordA[i] = mergeArray[i];

	delete [] mergeArray;

	//In final arrays shapes are blocks of 4 cordinates, 2 from each array, all shapes are independent and no shapes ceiling acts as anothers floor
	//A shapes floor cordinates are not actually included in the shape, so when we eventually pass to smerge, the floors will need a minus 1.
	//When a ceiling and a floor are the same, nothing is included, this will need to be another condition checked before passing to merge


	//----------MERGE SHAPES----------
	int batchNum2 = 0;
	//This makes think
	//Sends select values to different processors to do SRANK calculation
	int * windexPtr = &FinalWindex[0];
	while((my_rank+batchNum2) < int(numShapeCoord/2))
	{
		int i = 2*(my_rank+batchNum2);
		int winPos = FinalWindex[i/2]+1;
		if(shapeCordA[i]==shapeCordA[i+1])//Shape with and empty left side
		{
			int shapeSize = shapeCordB[i+1]-shapeCordB[i];
			int bValueIndex = shapeCordB[i];
			for(int x=0; x < shapeSize; x++)
			{
				localWin[winPos+x] = b[bValueIndex+x];
			}
		}
		if(shapeCordB[i]==shapeCordB[i+1])//Shape with and empty right side
		{
			int shapeSize = shapeCordA[i+1]-shapeCordA[i];
			int aValueIndex = shapeCordA[i];
			for(int x=0; x < shapeSize; x++)
			{
				localWin[winPos+x] = a[aValueIndex+x];
			}
		}
		int top1 = shapeCordA[i];
		int top2 = shapeCordB[i];
		int bottom1 = shapeCordA[i+1];
		int bottom2 = shapeCordB[i+1];
		//When shape sides actually need to be merged
		int * localOutput = new int[bottom1-top1+bottom2-top2];
		smerge(&a[top1],&b[top2],bottom1-top1-1,bottom2-top2-1, localOutput);
		for(int x=0; x < bottom1-top1+bottom2-top2 ; x++)
		{
			localWin[winPos+x] = localOutput[x];
		}

		batchNum2 += p;
	}

	MPI_Allreduce(localWin,&WIN[0],lasta+lastb+2,MPI_INT, MPI_MAX,MPI_COMM_WORLD);

	for(int i = 0; i < lasta+lastb+2; i++)
	{
		output[i] = WIN[i];
	}

	delete [] WIN;
	delete [] FinalWindex;
}
