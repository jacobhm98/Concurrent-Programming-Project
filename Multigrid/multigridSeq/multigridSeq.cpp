#include "../Jacobi.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

//function declarations
void initializeGrid(int gridSize, vector<vector<vector<double>>> &Matrix);
void resize(int gridSize, vector<vector<vector<double>>> &Matrix);
void setDirichletBoundaryConditions(vector<vector<vector<double>>> &Matrix);
void restrict(vector<vector<vector<double>>> &Matrix);

//global datastructure


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
	if ((gridSize - 5) % 4 != 0){
		cout << "The gridsize - 5 needs to be divisible by 4 for the reduction to work properly!" << endl;
	}
	
	vector<vector<vector<double>>> Matrix;
	initializeGrid(gridSize, Matrix);
	cout << "initialized" << endl;
	//begin the computations, start the timer right before
	auto startTime = std::chrono::high_resolution_clock::now();
	Jacobi iterate(Matrix);
	auto endTime = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
	//print out the specified data
	//printf("Grid Size, and Number of Iterations: %d, %d\n", gridSize, numIters);
	//cout << "Execution time of the computational part, in microseconds: " << duration.count() << endl;
	//cout << "The largest change an arbitrary grid went through this cycle is: " << maxDifference << endl;
	//print out the state of the matrix to filedata.out
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

//a method which takes the gridsize and initializes our matrix to what it's supposed to be (border cells = 1, everything else = 0)
void initializeGrid(int gridSize, vector<vector<vector<double>>> &Matrix){
	//initialize the global matrix
	resize(gridSize, Matrix);

	//set boundary values to 1
	setDirichletBoundaryConditions(Matrix);

	//set internal values to 0
	for (int i = 1; i < gridSize - 1; ++i){
		for (int j = 1; j < gridSize - 1; ++j){
			Matrix[0][i][j] = 0;
		}
	}
}

void restrict(vector<vector<vector<double>>> &Matrix){
	int newSize = (Matrix[0].size() + 1) / 2;
	vector<vector<vector<double>>> temp;
	

}
void resize(int gridSize, vector<vector<vector<double>>> &Matrix){
	Matrix.resize(2);
	Matrix[0].resize(gridSize);
	Matrix[1].resize(gridSize);
	for (int i = 0; i < gridSize; ++i){
		Matrix[0][i].resize(gridSize);
		Matrix[1][i].resize(gridSize);
	}
	
}
void setDirichletBoundaryConditions(vector<vector<vector<double>>> &Matrix){

	int gridSize = Matrix[0].size();
	for (int i = 0; i < gridSize; i += gridSize - 1){
		for (int j = 0; j < gridSize; ++j){
			Matrix[0][i][j] = 1;
			Matrix[0][j][i] = 1;
			Matrix[1][i][j] = 1;
			Matrix[1][j][i] = 1;
		} 
	}
}

