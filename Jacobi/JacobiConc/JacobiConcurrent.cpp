/* A sequential version of a PDE solver for Laplace's equation. Using the Jacobi iteration technique.
 * Compile: g++ JacobiConcurrent.cpp -o <executable>
 * Usage: ./<executable> int gridSize, int number of iterations, int number of workers
 */
#ifndef _REENTRANT
#define _REENTRANT
#endif
#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

//method declarations
void initializeGrid(int);
void updateMatrix(int, int);
double maxDiff(int);
void printMatrix(int, int);
void printMatrixtoFile(int current);

//global datastructure
vector<vector<vector<double>>> Matrix;


int main (int argc, char * argv[]){
	cout << std::fixed;
	cout << std::setprecision(2);
	int numIters = 0;
	int gridSize = 0;
	int numWorkers = 1;
	if (argc != 4){
		cout << "Arguments should be gridSize, numIters, and numWorkers" << endl;
		return 1;
	}
	else{
		gridSize =  atoi(argv[1]);
		numIters = atoi(argv[2]);
		numWorkers = atoi(argv[2]);
	}
	//initialize grid, set number of workers
	initializeGrid(gridSize);

	//begin the computations, start the timer right before
	auto startTime = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < numIters; ++i){
		updateMatrix(gridSize, numWorkers);
	}
	double maxDifference = maxDiff(gridSize);
	auto endTime = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime);
	//print out the specified data
	printf("Grid Size, and Number of Iterations: %d, %d\n", gridSize, numIters);
	cout << "Execution time of the computational part, in microseconds: " << duration.count() << endl;
	cout << "The maximum error is: " << maxDifference << endl;
	//print out the state of the matrix to filedata.out
	printMatrixtoFile(0);
	return 0;
}
//a method which takes the gridsize and initializes our matrix to what it's supposed to be (border cells = 1, everything else = 0)
void initializeGrid(int gridSize){
	//initialize the global matrix
	Matrix.resize(2);
	Matrix[0].resize(gridSize);
	Matrix[1].resize(gridSize);
	for (int i = 0; i < gridSize; ++i){
		Matrix[0][i].resize(gridSize);
		Matrix[1][i].resize(gridSize);
	}

	//set boundary values to 1
	for (int i = 0; i < gridSize; i += gridSize - 1){
		for (int j = 0; j < gridSize; ++j){
			Matrix[0][i][j] = 1;
			Matrix[0][j][i] = 1;
			Matrix[1][i][j] = 1;
			Matrix[1][j][i] = 1;
		} 
	}
	//set internal values to 0
	for (int i = 1; i < gridSize - 1; ++i){
		for (int j = 1; j < gridSize - 1; ++j){
			Matrix[0][i][j] = 0;
		}
	}
}
//Laplace's PDE, or Nabla² takes the limit of each point tending towards the average of its neighbours points, 
//and how this varies with time. This is why it is sometimes called the heat equation: given an input in 1d, 2d, 3d ... 
//where each point has a value corresponding to its temperature, where the temperature is given as f(x,y,z..) , as time -> infinity,
//the temperatures will reach an equilibrium across the entire input space. A hot point surrounded by very cold points will get colder quickly
//and the cold points will get slightly warmer, until all of the points reach the same value. Laplaces PDE: Nabla²(f(x,y,z...) = 0.
//This method simulates this, by updating each point on the matrix to the average of its neighbouring 4 points.
void updateMatrix(int gridSize, int numWorkers){
	//set the number of threads we wish to use, taken from the commandline. Do the parallelize the for loop with omp.	
	omp_set_num_threads(numWorkers);
	#pragma omp parallel for
	//do the calculations for all internal points of the grid (not boundary points)
	for (int i = 1; i < gridSize - 1; ++i){
		for (int j = 1; j < gridSize - 1; ++j){
			Matrix[1][i][j] = (Matrix[0][i-1][j]+Matrix[0][i+1][j]+Matrix[0][i][j-1]+Matrix[0][i][j+1]) * .25;
		}
	}
	
	//make new matrix old matrix to prepare for next iteration, save old matrix in new matrix so we can calculate maxdiff still
	vector<vector<double>> temp;
	vector<vector<double>> * p1 = &Matrix[0];
	vector<vector<double>> * p2 = &Matrix[1];
	temp = *p2;
	*p2 = *p1;
	*p1 = temp;
}
//Calculate the values delta of each point from this iteration and last, and return the value of the highest one. The smaller this value is,
//the less effective each subsequent iteration is, and the closer we are to the limit.
double maxDiff(int gridSize){
	double maxDiff = 0;
	for (int i = 0; i < gridSize; ++i){
		for (int j = 0; j < gridSize; ++j){
			double currentDiff = (1 - Matrix[0][i][j]);
			if (maxDiff < currentDiff){
				maxDiff = currentDiff;
			}
			}
		}
	return maxDiff;
}

//For debugging.
void printMatrix(int gridSize, int current){
	for (int i = 0; i < gridSize; ++i){
		for (int j = 0; j < gridSize; ++j){
			cout << "|";
			cout << Matrix[current][i][j];
			cout << "|";
		}
		cout << endl;
	} 
}
void printMatrixtoFile(int current){
	std::ofstream out;
	out << std::fixed << std::setprecision(2);
	out.open("./filedata.out");
	int gridSize = Matrix[0].size();
	for (int i = 0; i < gridSize; ++i){
		for (int j = 0; j < gridSize; ++j){
			out << "|";
			out << Matrix[current][i][j];
			out << "|";
		}
		out << endl;
	} 
}
