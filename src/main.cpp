#include "main.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "point.h"
#include "io.h"
#include "maths.h"
#include "params.h"
#include "sausages.h"

using namespace std;
using namespace model;

// Declare model variables (need to be extern in main.h)
///@TODO Global variables are bad, make this a class
std::vector<vector3d> model::colloidPos; ///< Positions of the colloids
std::vector<double> model::colloidRadii; ///< Radii of the colloids
std::vector<Point> model::allPoints; ///< All points in simulation
std::vector<Sausage> model::allSausages; ///< All sausages in simulation

/**
 * Main function
 * Program takes input filename as first argument
 * Optional second argument is a JSON param file
 */
int main(int argc, char *argv[]){

	// Check #args
	if (argc!=2 && argc!=3){
		cerr<<"Usage: "<<argv[0]<<" inputFile [paramFile]"<<endl;
		exit(EXIT_FAILURE);
	}

	// Set various parameters and variables
	if (argc>2){
		set_params(argv[2]);
	} else {
		// Just use defaults
		//set_params(NULL);
	}
	debug() << "Colloids in param file: " << params::colloidsInParamFile << endl;

	// Read input data file
	info() << "Reading file " << argv[1] << endl;
	int num_below_threshold = read_input(argv[1],allPoints);
	info() << "Read "<<num_below_threshold << " points."<<endl;
	info() << "Stored "<< allPoints.size() << " points."<<endl;

	// Separate the points into separate, contiguous sausages
	info() << "Distinguishing sausages..." << endl;
	flood_fill_separate(allPoints,allSausages);

	verbose() << "Sausage sizes:" << endl;
	for (size_t i=0; i<allSausages.size(); i++){
		verbose() << "  " << i << " " << allSausages[i].points.size() << endl;
	}

	// If the sausages are 'too small', ignore them.
	vector<int> relevant_sausages; ///< A vector of sausageIDs of 'sufficiently large' sausages
	for (size_t i=0; i<allSausages.size(); i++){
		int sausageSize=allSausages[i].points.size();
		if (sausageSize < params::silent_ignore_size*num_below_threshold ){
			verbose() << "Ignoring sausage #" << i << endl;
			allSausages[i].is_significant=false;
		} else if (sausageSize > params::silent_ignore_size*num_below_threshold &&
			sausageSize < params::min_sausage_size*num_below_threshold   ){
			warning() << "Sausage #" << i << " of size "  <<
				sausageSize << " is not tiny, but is being ignored." << endl;
			allSausages[i].is_significant=false;
		} else if (sausageSize > params::min_sausage_size*num_below_threshold){
			allSausages[i].is_significant=true;
			relevant_sausages.push_back(i);
		}
	}
	info() << "Found " << relevant_sausages.size() << " sufficiently large sausages." << endl;

	info() << "Finding endpoints..."<<endl;
	for (size_t i=0;i<relevant_sausages.size();i++){
		allSausages[relevant_sausages[i]].find_endpoints();
	}

	info() << "Joining small gaps..."<<endl;
	join_endpoints(allSausages,relevant_sausages);
	debug() << "Sausage sizes:" << endl;
	for (size_t i=0; i<allSausages.size(); i++){
		debug() << "  " << i << " " << allSausages[i].points.size() << endl;
	}
	debug() << "Relevant sausages: " << endl;
	for (size_t i=0; i<relevant_sausages.size(); i++){
		debug() << "  " << i << " " << relevant_sausages[i] << endl;
	}

	// Exit program is number of sausages found is not equal 2 or 1.
	if ( relevant_sausages.size() > 2) {
		cerr<<"Number of sausage sizes detected is larger than 2! Unphysical. Check your input parameters." << endl;
		exit(EXIT_FAILURE);
	}

	// Analysis for two relevant sausages found 
	else if ( relevant_sausages.size() == 2) {
	double size [2];
	double ratio;
		info() << "Two relevant sausages found." << endl;
		debug()<<"#points: "<<allSausages[relevant_sausages[0]].points.size()<<" and "<<allSausages[relevant_sausages[1]].points.size()<<endl;
		for (size_t i=0; i<relevant_sausages.size(); i++){
			size[i] = allSausages[relevant_sausages[i]].points.size();
			allSausages[relevant_sausages[i]].find_com();
			allSausages[relevant_sausages[i]].find_pobf();
			allSausages[relevant_sausages[i]].estimate_sausage_length();
			info() << "Size of i th sausage: " << i << " " << size[i] << endl;
		}
		ratio = fabs(2*(size[0]-size[1])/(size[0]+size[1]));
		if (ratio < params::ratio_two_rings ){
			info() << "The sausages have very similar size. Hence it is the Two-ring structure." << endl;
		} else if (ratio > params::ratio_2nd_loop ) {
			info() << "The sausages have very different size. Hence it is the 2nd_loop structure." << endl;
		} else {
			cerr << "Two relevant sausages were found but the ratio was not distinctive enought to distinguish between Two-ring and 2nd_loop. \
			Maybe change input parameters." << endl;
			exit(EXIT_FAILURE);
			}
	}
	
	// Analysis for one relevant sausages found 
	else if ( relevant_sausages.size() == 1) {
		info() << "One relevant sausage found." << endl;
		double size = allSausages[relevant_sausages[0]].points.size();
		allSausages[relevant_sausages[0]].find_com();
		allSausages[relevant_sausages[0]].find_pobf();
		//allSausages[relevant_sausages[0]].track_sausage();
		info() << "Size of sausage: " << size << endl;
		//Do FF to distinguish between figure of eight and figure of omega 
		allSausages[relevant_sausages[0]].flood_fill_classify();
	}

	// Wrap up and exit
	info() << "Exiting successfully" << endl;
	exit(EXIT_SUCCESS);
}


