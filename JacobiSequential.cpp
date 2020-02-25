#include <iostream>
#include <cstdlib>
#include <fstream>
#include <chrono>
#define DEBUG 0

using std::cout;
using std::endl;

//method declarations
void initializeGrid(int);
void updateMatrix(int);
int maxDiff(int);

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
	auto startTime = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < numIters; ++i){
		updateMatrix(gridSize);
	}
	int maxDifference = maxDiff(gridSize);
	auto endTime = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
	printf("Grid Size, and Number of Iterations: %d, %d\n", gridSize, numIters);
	cout << "Execution time of the computational part, in microseconds: " << duration.count() << endl;
	cout << "The largest change an arbitrary grid went through this cycle is: " << maxDifference << endl;
	std::ofstream out;
	out.open("./filedata.out");
	for (int i = 0; i < gridSize; ++i){
		for (int j = 0; j < gridSize; ++j){
			out << "|" << Matrix[0][i][j] << "|";
		}
		out << endl;
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

int maxDiff(int gridSize){
	double maxDiff = 0;
	for (int i = 0; i < gridSize; ++i){
		for (int j = 0; j < gridSize; ++j){
			int currentDiff = std::abs(Matrix[0][i][j] - Matrix[1][i][j]);
			if (maxDiff < currentDiff){
				maxDiff = currentDiff;
			}
		}
	}
	return maxDiff;
}
