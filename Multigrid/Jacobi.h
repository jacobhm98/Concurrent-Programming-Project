
#pragma once
#include <vector>

using std::vector;

class Jacobi {
	void updateMatrix(vector<vector<vector<double>>> &Matrix);
	double maxDiff(vector<vector<vector<double>>> &Matrix);
	public:
	void iterate(int numIters, vector<vector<vector<double>>> &Matrix);
	void iterate(vector<vector<vector<double>>> &Matrix);
	void printMatrix(vector<vector<vector<double>>> &Matrix, int);
};
