#include "../Jacobi.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <vector>

using std::cout;
using std::endl;
using std::vector;
//global datastructure
vector<vector<vector<double>>> Matrix;


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
	
	initializeGrid(gridSize);
	//begin the computations, start the timer right before
	auto startTime = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < numIters; ++i){
		updateMatrix(gridSize);
	}
	double maxDifference = maxDiff(gridSize);
	auto endTime = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
	//print out the specified data
	printf("Grid Size, and Number of Iterations: %d, %d\n", gridSize, numIters);
	cout << "Execution time of the computational part, in microseconds: " << duration.count() << endl;
	cout << "The largest change an arbitrary grid went through this cycle is: " << maxDifference << endl;
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
