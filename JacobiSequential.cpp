#include <iostream>
#define DEBUG 0

using std::cout;
using std::endl;

//method declarations
void initializeGrid(int);
void updateMatrix(int);

//global datastructures
int*** Matrix;


int main (int argc, char * argv[]){
	int numIters = 0;
	int gridSize = 0;
	if (argc != 3){
		cout << "Arguments should be gridSize, and numIters" << endl;
		return 1;
	}
	else{
		gridSize =  atoi(argv[1]);
		numIters = atoi(argv[2]);
	}
	for (int i = 0; i < numIters; ++i){
		updateMatrix(gridSize);
	}	

	return 0;
}

void initializeGrid(int gridSize){
	//initialize the global matrix
	//Matrix[0] is old matrix, Matrix[1] is new matrix.
	Matrix = new int**[2];
	Matrix[0] = new int*[gridSize];
	Matrix[1] = new int*[gridSize];
	for (int i = 0; i<gridSize; ++i){
		Matrix[0][i], Matrix[1][i] = new int[gridSize];
	}

	//set boundary values to 1
	for (int i = 0; i < gridSize; i += gridSize - 1){
		for (int j = 0; j < gridSize; ++j){
			Matrix[0][i][j] = 1;
			Matrix[0][j][i] = 1;
		} 
	}
	//set internal values to 0
	for (int i = 1; i < gridSize - 1; ++i){
		for (int j = 1; j < gridSize - 1; ++j){
			Matrix[0][i][j] = 0;
		}
	}
}

void updateMatrix(int gridSize){
	//do the calculations for all internal points of the grid (not boundary points)
	for (int i = 1; i < gridSize - 1; ++i){
		for (int j = 1; j < gridSize - 1; ++j){
			Matrix[1][i][j] = (Matrix[0][i-1][j]+Matrix[0][i+1][j]+Matrix[0][i][j-1]+Matrix[0][i][j+1]) * .25;
		}
	}
	//make new matrix old matrix to prepare for next iteration, save old matrix in new matrix so we can calculate maxdiff still
	int ** temp = Matrix[0];
	Matrix[0] = Matrix[1];
	Matrix[1] = temp;
}
