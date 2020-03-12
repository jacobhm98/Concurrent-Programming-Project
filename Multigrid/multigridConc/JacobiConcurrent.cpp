
#include <omp.h>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "../Jacobi.h"
using std::cout;
using std::endl;
using std::vector;

//iterate number passed in for the coarse grid
void Jacobi::iterateP(int numWorkers, int numIters, vector<vector<vector<double>>> &Matrix){
	for (int i = 0; i < numIters; ++i){
		updateMatrix(numWorkers, Matrix);
	}
}

//iterate four times on the finer grid
void Jacobi::iterateP(int numWorkers, vector<vector<vector<double>>> &Matrix){
	for (int i = 0; i < 4; ++i){
		updateMatrix(numWorkers, Matrix);
	}
}

//Laplace's PDE, or Nabla² takes the limit of each point tending towards the average of its neighbours points, 
//and how this varies with time. This is why it is sometimes called the heat equation: given an input in 1d, 2d, 3d ... 
//where each point has a value corresponding to its temperature, where the temperature is given as f(x,y,z..) , as time -> infinity,
//the temperatures will reach an equilibrium across the entire input space. A hot point surrounded by very cold points will get colder quickly
//and the cold points will get slightly warmer, until all of the points reach the same value. Laplaces PDE: Nabla²(f(x,y,z...) = 0.
//This method simulates this, by updating each point on the matrix to the average of its neighbouring 4 points.
void Jacobi::updateMatrix(int numWorkers, vector<vector<vector<double>>> &Matrix){
	int gridSize = Matrix[0].size();
	
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

//calculate the difference between 1 and the smallest valued grid point
double Jacobi::maxDiff(vector<vector<vector<double>>> &Matrix){
	int gridSize = Matrix[0].size();
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
void Jacobi::printMatrix(vector<vector<vector<double>>> &Matrix, int current){
	cout << "hello" << endl;
	int gridSize = Matrix[0].size();
	for (int i = 0; i < gridSize; ++i){
		for (int j = 0; j < gridSize; ++j){
			cout << "|";
			cout << Matrix[current][i][j];
			cout << "|";
		}
		cout << endl;
	} 
}
