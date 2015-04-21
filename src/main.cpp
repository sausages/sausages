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

	// Setup model
	Model model;

	// Set various parameters and variables
	if (argc>2){
		set_params(argv[2],model.colloidPos);
	} else {
		// Just use defaults
		set_params(NULL,model.colloidPos);
	}
	brief({1}) << "brief_version	1" << endl;
	brief({2}) << "brief_version	2" << endl;
	brief({1,2}) << "input_filename	" << argv[1] << endl;

	debug() << "Colloids in param file: " << params::colloidsInParamFile << endl;

	// Read input data file
	info() << "Reading file " << argv[1] << endl;
	read_input(argv[1], model);
	info() << "Read "<< model.total_points << " points."<<endl;
	info() << "Stored "<< model.allPoints.size() << " points."<<endl;

	brief({1,2}) << "total_points	" << model.total_points << endl;
	brief({1,2}) << "threshold	" << params::threshold << endl;
	brief({1,2}) << "threshold_pts	" << model.threshold_points.size() << endl;

	if (!params::thresholded_filename.empty()) write_points(params::thresholded_filename, model.allPoints, model.threshold_points);

	// Separate the points into separate, contiguous sausages
	info() << "Distinguishing sausages..." << endl;
	flood_fill_separate(model.allPoints,model.allSausages);

	verbose() << "Sausage sizes:" << endl << flush;
	for (size_t i=0; i<model.allSausages.size(); i++){
		verbose() << "  " << i << " " << model.allSausages[i].points.size() << endl;
	}

	if (!params::prejoin_sausage_filename.empty()){
		debug() << "Printing pre-joined sausages to files: " << params::prejoin_sausage_filename << endl;
		for (size_t i=0; i<model.allSausages.size(); i++){
			write_points(params::prejoin_sausage_filename+to_string(i), model.allSausages[i].points);
		}
	}

	// Find endpoints of all relevant sausages
	info() << "Finding endpoints..."<<endl;
	for (size_t i=0;i<model.allSausages.size();i++){
		model.allSausages[i].find_endpoints();
	}

	// Join the endpoints
	debug() << "We currently have " << model.allSausages.size() << " sausages" << endl;
	info() << "Joining small gaps..."<<endl;
	//join_endpoints(model.allSausages,relevant_sausages);
	join_endpoints(model);
	debug() << "Now we have " << model.allSausages.size() << " sausages" << endl;


	// If the sausages are 'too small', ignore them.
	vector<int> relevant_sausages; ///< A vector of sausageIDs of 'sufficiently large' sausages
	for (size_t i=0; i<model.allSausages.size(); i++){
		int sausageSize=model.allSausages[i].points.size();
		if (sausageSize < params::silent_ignore_size*model.threshold_points.size() ){
			verbose() << "Ignoring sausage #" << i << endl;
			model.allSausages[i].is_significant=false;
		} else if (sausageSize > params::silent_ignore_size*model.threshold_points.size() &&
			sausageSize < params::min_sausage_size*model.threshold_points.size()   ){
			warning() << "Sausage #" << i << " of size "  <<
			sausageSize << " is not tiny, but is being ignored." << endl;
			model.allSausages[i].is_significant=false;
		} else if (sausageSize > params::min_sausage_size*model.threshold_points.size()){
			model.allSausages[i].is_significant=true;
			relevant_sausages.push_back(i);
		}
	}
	info() << "Found " << relevant_sausages.size() << " sufficiently large sausages." << endl;

	if (!params::sausage_filename.empty()){
		debug() << "Printing joined sausages to files: " << params::sausage_filename << endl;
		for (size_t i=0; i<model.allSausages.size(); i++){
			write_points(params::sausage_filename+to_string(i), model.allSausages[i].points);
		}
	}

	// Check that all relevant sausages are closed loops, open loops are unphysical, except along pbc
	info() << "Checking all relevant sausages are closed loops..." << endl;
	for (size_t i=0; i<relevant_sausages.size(); i++){
        model.allSausages[relevant_sausages[i]].flood_fill_closed_loops(model.colloidPos);
    }
	info() << "Check successful." << endl;

	// Plenty of output
	brief({1,2}) << "num_sausages	" << model.allSausages.size() << endl;

	debug() << "Post-joining sausage sizes:" << endl;
	for (size_t i=0; i<model.allSausages.size(); i++){
		debug() << "  " << i << " " << model.allSausages[i].points.size() << endl;
		brief({1,2}) << "saus_" << i << "_size	" << model.allSausages[i].points.size() << endl;
	}
	debug() << "Relevant sausages (prints index of sausages): " << endl;
	for (size_t i=0; i<relevant_sausages.size(); i++){
		debug() << "  " << i << " " << relevant_sausages[i] << endl;
	}

	// Exit program is number of sausages found is not equal 2 or 1.
	if ( relevant_sausages.size() > 2) {
		cerr<<"More than two sausages found after endpoint-joining! This is unphysical, please check your input parameters." << endl;
		exit(EXIT_FAILURE);
	}

	// Analysis for two relevant sausages found 
	else if ( relevant_sausages.size() == 2) {
		double size [2];
		double ratio;
		info() << "Two relevant sausages found." << endl;
		debug()<<"#points: "<<model.allSausages[relevant_sausages[0]].points.size()<<" and "<<model.allSausages[relevant_sausages[1]].points.size()<<endl;
		for (size_t i=0; i<relevant_sausages.size(); i++){
			size[i] = model.allSausages[relevant_sausages[i]].points.size();
			model.allSausages[relevant_sausages[i]].find_com();
			model.allSausages[relevant_sausages[i]].find_pobf();
			model.allSausages[relevant_sausages[i]].sphere_tracking();
            model.allSausages[relevant_sausages[i]].calculate_sausage_length_spheres();
			info() << "Size of i th sausage: " << i << " " << size[i] << endl;
		}
	// Find the ratio in size for the two sausages. 
		ratio = fabs(2*(size[0]-size[1])/(size[0]+size[1]));
		if (ratio < params::ratio_two_rings ){
			info() << "The sausages have very similar size. Hence it is the Two-ring structure." << endl;
			brief({1,2}) << "system_class	3" << endl;
			brief({1,2}) << "system_str	twoRings" << endl;
		} else if (ratio > params::ratio_2nd_loop ) {
			info() << "The sausages have very different size. Hence it is the 2nd_loop structure." << endl;
			brief({1,2}) << "system_class	4" << endl;
			brief({1,2}) << "system_str	theta" << endl;
		} else {
			cerr << "Two relevant sausages were found but the ratio was not distinctive enought to distinguish between Two-ring and 2nd_loop. \
			Maybe change input parameters." << endl;
			brief({1,2}) << "system_class	0" << endl;
			brief({1,2}) << "system_str	undefined" << endl;
			exit(EXIT_FAILURE);
			}
	}

	// Analysis for one relevant sausages found 
	else if ( relevant_sausages.size() == 1) {
		info() << "One relevant sausage found." << endl;
		double size = model.allSausages[relevant_sausages[0]].points.size();
		model.allSausages[relevant_sausages[0]].find_com();
		model.allSausages[relevant_sausages[0]].find_pobf();
		model.allSausages[relevant_sausages[0]].sphere_tracking();
		model.allSausages[relevant_sausages[0]].calculate_sausage_length_spheres();
		info() << "Size of sausage: " << size << endl;
		//Do FF to distinguish between figure of eight and figure of omega 
		int classification = model.allSausages[relevant_sausages[0]].flood_fill_classify(model.colloidPos);
		brief({1}) << "system_class	" << classification << endl;
		switch (classification){
			case 1:	brief({1}) << "system_str	figureEightLHS" << endl; break;
			case 2:	brief({1}) << "system_str	figureEightRHS" << endl; break;
			case 5:	brief({1}) << "system_str	omega" << endl; break;
			default: error() << "Unexpected system classification" << endl; exit(EXIT_FAILURE);
		}
	}

	// Wrap up and exit
	for (Point* p : model.allPoints) {delete p;}
	model.allSausages.clear();
	model.allPoints.clear();

	info() << "Exiting successfully" << endl;
	return EXIT_SUCCESS;
}

