#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "point.h"
#include "io.h"
#include "find_sausages.h"
#include "main.h"

using namespace std;
using namespace params;

verbosityLevel params::verbosity=NORMAL;
const double params::threshold_level = 0.12;

/**
 * Main function
 * Program takes input filename as first (and only) argument
 */
int main(int argc, char *argv[]){
	// Set various parameters and variables
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

	if (verbosity>=NORMAL) cout << "Reading file " << argv[1] << endl;
	read_xyzclcpcs(infile,allPoints);


	// Put all points with cl<threshold in a sausage
	if (verbosity>=NORMAL) cout << "Thresholding, sausages have cl<" << params::threshold_level << endl;
	threshold(allPoints);

	// Count the sausages
	if (verbosity>=NORMAL) cout << "Pixel counts in (0) or not in (1) a sausage:" << endl;
	vector<int> sausage_count=count_sausages(allPoints);
	for (size_t i=0; i<sausage_count.size(); i++){
		cout << i << " " << sausage_count[i] << endl;
	}

	// Link the points to their neighbours
	if (verbosity>=NORMAL) cout << "Linking pixels..." << endl;
	neighLink_xyzclcpcs(allPoints);

	// Check the neighbours
	if (verbosity>=DEBUG) printAllNeighs(allPoints);

	// Separate the points into separate, contiguous sausages
	if (verbosity>=NORMAL) cout << "Distinguishing sausages..." << endl;
	flood_fill(allPoints);

	// Count the sausages
	sausage_count=count_sausages(allPoints);
	if (verbosity>=NORMAL){
		cout << "Sausage sizes:" << endl;
		for (size_t i=0; i<sausage_count.size(); i++){
			cout << i << " " << sausage_count[i] << endl;
		}
	}

	if (verbosity>=NORMAL) cout << "Exiting successfully" << endl;
}

