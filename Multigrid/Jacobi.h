#pragma once

class Jacobi {

	public:
	//constructors
	Jacobi(int numIters, vector<vector<vector<double>>> &Matrix);
	Jacobi(vector<vector<vector<double>>> &Matrix);
	//method declarations
	void initializeGrid(int);
	void updateMatrix(int);
	double maxDiff(int);
	void printMatrix(int, int);
}
