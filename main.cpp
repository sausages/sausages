#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "main.h"
#include "io.h"

using namespace std;

/**
 * Main function
 * Program takes input filename as first (and only) argument
 */
int main(int argc, char *argv[]){
	if (argc!=2){
		cerr<<"Usage: "<<argv[0]<<" inputFile"<<endl;
		exit(EXIT_FAILURE);
	}
	ifstream infile (argv[1]);
	if (!infile.is_open()){
		cerr<<"Error opening file "<<argv[1]<<endl;
		exit(EXIT_FAILURE);
	}

	vector<Point> allPoints; ///< A vector of all points in simulation

	read_xyzclcpcs(infile,allPoints);

	cout << allPoints[1].z << endl;

	cout << "Exiting successfully" << endl;
}




