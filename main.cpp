#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "point.h"
#include "io.h"
#include "find_sausages.h"

using namespace std;

/**
 * Main function
 * Program takes input filename as first (and only) argument
 */
int main(int argc, char *argv[]){
	// Set various parameters and variables
	const double threshold_level = 0.12;
	vector<Point> allPoints; ///< A vector of all points in simulation

	// Check #args
	if (argc!=2){
		cerr<<"Usage: "<<argv[0]<<" inputFile"<<endl;
		exit(EXIT_FAILURE);
	}

	// Open and read input data file
	ifstream infile (argv[1]);
	if (!infile.is_open()){
		cerr<<"Error opening file "<<argv[1]<<endl;
		exit(EXIT_FAILURE);
	}
	read_xyzclcpcs(infile,allPoints);


	// Put all points with cl>threshold in a sausage
	threshold(allPoints,threshold_level);

	// Count the sausages
	vector<int> sausage_count=count_sausages(allPoints);
	for (vector<int>::const_iterator it=sausage_count.begin();it!=sausage_count.end(); ++it){
		cout << *it << endl;
	}

	cout << "Exiting successfully" << endl;
}

