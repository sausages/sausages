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

namespace params{
	const verbosityLevel verbosity = NORMAL;
	const double threshold_level = 0.12;
	const double silent_ignore_size = 0.01;
	const double min_sausage_size = 0.1;
}

/**
 * Main function
 * Program takes input filename as first (and only) argument
 */
int main(int argc, char *argv[]){
	// Set various parameters and variables
	vector<Point> allPoints; ///< A vector of all points in simulation
	vector<int> sausage_count; ///< A vector of the size of all sausages
	vector<int> relevant_sausages; ///< A vector of sausageIDs of 'sufficiently large' sausages

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
	if (verbosity>=NORMAL) cout << "  found " << allPoints.size() << " points." << endl;

	// Put all points with cl<threshold in a sausage
	if (verbosity>=NORMAL) cout << "Thresholding, sausages have cl<" << params::threshold_level << endl;
	int num_below_threshold = threshold(allPoints);
	if (verbosity>=NORMAL) cout << "  " << num_below_threshold << " out of "
		<< allPoints.size() << " points were below the theshold." << endl;

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
			cout << "  " << i << " " << sausage_count[i] << endl;
		}
	}

	// Measure the sausages' sizes. If they are 'too small', ignore them.
	for (size_t i=2; i<sausage_count.size(); i++){
		if (sausage_count[i] < params::silent_ignore_size*num_below_threshold &&
			verbosity >= VERBOSE ){
			cout << "Ignoring sausage #" << i << endl;
		} else if (sausage_count[i] > params::silent_ignore_size*num_below_threshold &&
			sausage_count[i] < params::min_sausage_size*num_below_threshold   ){
			if (verbosity>=WARNING) cout << "Sausage #" << i << " of size "  <<
				sausage_count[i] << " is not tiny, but is being ignored." << endl;
		} else if (sausage_count[i] > params::min_sausage_size*num_below_threshold){
			relevant_sausages.push_back(i);
		}
	}
	if (verbosity>=NORMAL)
		cout << "Found " << relevant_sausages.size() << " sufficiently large sausages." << endl;

	// Wrap up and exit
	if (verbosity>=NORMAL) cout << "Exiting successfully" << endl;
}

