
#ifndef _REENTRANT
#define _REENTRANT
#endif
#include <omp.h>
#include <iostream>
#define DEBUG 0
#define DEFAULT_THREADS 1

using std::cout;
using std::endl;

int main (int argc, char * argv[]){
	if (argc != 4){
		cout << "Arguments should be gridSize, numIters, and numWorkers" << endl;
		return 1;
	}
	else{
		int gridSize =  atoi(argv[1]);
		int numIters = atoi(argv[2]);
		int numWorkers = atoi(argv[3]);
		cout << "test" << endl;
	}
	return 0;
}
