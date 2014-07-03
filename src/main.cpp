#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "point.h"
#include "io.h"
#include "sausages.h"
#include "main.h"

using namespace std;

namespace params{
	//const verbosityLevel verbosity = INFO;
	const verbosityLevel verbosity = INFO;
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
	vector<Sausage> allSausages; ///< A vector of all sausages in simulation
	vector<int> relevant_sausages; ///< A vector of sausageIDs of 'sufficiently large' sausages

	// Check #args
	if (argc!=2){
		cerr<<"Usage: "<<argv[0]<<" inputFile"<<endl;
		exit(EXIT_FAILURE);
	}

	// Open, read and close input data file
	ifstream infile (argv[1]);
	if (!infile.is_open()){
		cerr<<"Error opening file "<<argv[1]<<endl;
		exit(EXIT_FAILURE);
	}

	info() << "Reading file " << argv[1] << endl;
	read_xyzclcpcs(infile,allPoints);
	info() << "  found " << allPoints.size() << " points." << endl;
	infile.close();

	// Put all points with cl<threshold in a sausage
	info() << "Thresholding, sausages have cl<" << params::threshold_level << endl;
	int num_below_threshold = threshold(allPoints);
	info() << "  " << num_below_threshold << " out of "
		<< allPoints.size() << " points were below the theshold." << endl;

	// Link the points to their neighbours
	info() << "Linking pixels..." << endl;
	neighLink_xyzclcpcs(allPoints);

	// Check the neighbours
	if (params::verbosity >= DEBUG) printAllNeighs(allPoints);

	// Separate the points into separate, contiguous sausages
	info() << "Distinguishing sausages..." << endl;
	flood_fill_separate(allPoints,allSausages);
	info() << "Sausage sizes:" << endl;
	for (size_t i=0; i<allSausages.size(); i++){
		info() << "  " << i+2 << " " << allSausages[i].points.size() << endl;
		}

	// If the sausages are 'too small', ignore them.
	for (size_t i=0; i<allSausages.size(); i++){
		int sausageSize=allSausages[i].points.size();
		if (sausageSize < params::silent_ignore_size*num_below_threshold ){
			verbose() << "Ignoring sausage #" << i+2 << endl;
			allSausages[i].is_significant=false;
		} else if (sausageSize > params::silent_ignore_size*num_below_threshold &&
			sausageSize < params::min_sausage_size*num_below_threshold   ){
			warning() << "Sausage #" << i+2 << " of size "  <<
				sausageSize << " is not tiny, but is being ignored." << endl;
			allSausages[i].is_significant=false;
		} else if (sausageSize > params::min_sausage_size*num_below_threshold){
			allSausages[i].is_significant=true;
			relevant_sausages.push_back(i);
		}
	}
	info() << "Found " << relevant_sausages.size() << " sufficiently large sausages." << endl;

	// Calculate properties of relevant sausages
	for (size_t i=0; i<relevant_sausages.size(); i++){
		Sausage thisSausage=allSausages[relevant_sausages[i]];
		thisSausage.find_com();
		thisSausage.find_pobf();
	}

	// Find centre of masses and orientation vectors following sausage
	for (size_t i=0; i<allSausages.size(); i++){
		if ( allSausages[i].is_significant==true ){
		allSausages[i].estimate_sausage_length();}
	}


	// Wrap up and exit
	info() << "Exiting successfully" << endl;
}


