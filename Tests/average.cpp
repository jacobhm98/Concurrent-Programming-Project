#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char ** argv){
	string file = argv[1];	
	ifstream in;
	ofstream out;
	in.open(file);
	int average = 0;
	string number;
	while (in >> number){
		average += stoi(number);
	}
	out.open(file);
	out << average/5;
	cout << average/5;

}
