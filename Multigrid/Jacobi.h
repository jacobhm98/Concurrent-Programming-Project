
#pragma once
#include <vector>

using std::vector;

class Jacobi {
	void updateMatrix(vector<vector<vector<double>>> &Matrix);
	double maxDiff(vector<vector<vector<double>>> &Matrix);
	void printMatrix(vector<vector<vector<double>>> &Matrix, int);
	public:
	//constructors
	Jacobi(int numIters, vector<vector<vector<double>>> &Matrix);
	Jacobi(vector<vector<vector<double>>> &Matrix);
	//method declarations
};
