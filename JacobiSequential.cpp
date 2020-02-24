
#ifndef _REENTRANT
#define _REENTRANT
#endif
#include <omp.h>
#include <iostream>
#define DEBUG 0
#define DEFAULT_THREADS 1

using std::cout;
using std::endl;

//method declarations
void initializeGrid(int);

//global datastructures
int** oldMatrix;
int** newMatrix;
int main (int argc, char * argv[]){
	if (argc != 4){
		cout << "Arguments should be gridSize, numIters, and numWorkers" << endl;
		return 1;
	}
	else{
		int gridSize =  atoi(argv[1]);
		int numIters = atoi(argv[2]);
		int numWorkers = atoi(argv[3]);
	}
	

	return 0;
}

void initializeGrid(int gridSize){
	//initialize the global matrix
	oldMatrix, newMatrix = new int*[gridSize];
	for (int i = 0; i<gridSize; ++i){
		oldMatrix[i], newMatrix[i] = new int[gridSize];
	}

	//set boundary values to 1
	for (int i = 0; i < gridSize; i += gridSize - 1){
		for (int j = 0; j < gridSize; ++j){
			oldMatrix[i][j] = 1;
			oldMatrix[j][i] = 1;
		} 
	}
	//set internal values to 0
	for (int i = 1; i < gridSize - 1; ++i){
		for (int j = 1; j < gridSize - 1; ++j){
			oldMatrix[i][j] = 0;
		}
	}
}
