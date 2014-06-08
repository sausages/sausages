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

	cout << "Reading file " << argv[1] << endl;
	read_xyzclcpcs(infile,allPoints);


	// Put all points with cl<threshold in a sausage
	cout << "Thresholding, sausages have cl<" << threshold_level << endl;
	threshold(allPoints,threshold_level);

	// Count the sausages
	cout << "Pixel counts in (0) or not in (1) a sausage:" << endl;
	vector<int> sausage_count=count_sausages(allPoints);
	for (size_t i=0; i<sausage_count.size(); i++){
		cout << i << " " << sausage_count[i] << endl;
	}

	// Link the points to their neighbours
	cout << "Linking pixels..." << endl;
	neighLink_xyzclcpcs(allPoints);

	// Check the neighbours
	//printAllNeighs(allPoints);

	// Separate the points into separate, contiguous sausages
	cout << "Distinguishing sausages..." << endl;
	flood_fill(allPoints);

	// Count the sausages
	cout << "Sausage sizes:" << endl;
	sausage_count=count_sausages(allPoints);
	for (size_t i=0; i<sausage_count.size(); i++){
		cout << i << " " << sausage_count[i] << endl;
		cout << endl;
	}

	cout << "Exiting successfully" << endl;
}

