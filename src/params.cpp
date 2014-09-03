#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include "cJSON/cJSON.h"
#include "io.h"
#include "params.h"

using namespace params;

// Set defaults
namespace params{
	verbosityLevel verbosity = INFO;
	double threshold_high = 0.12; // Sausages should be without gaps at this level
	double threshold_low = 0.04; // All different defect structures should be distinguishable (i.e. no ambigious blobs)
	double silent_ignore_size = 0.01; // If a sausage is smaller than this fraction of all points below the threshold, silently ignore it
	double min_sausage_size = 0.1; // A sausage is only 'significant' if it is larger than this fraction of all points below the threshold
	double pixel_size = 0.5; // Distance between nearest-neighbour points
	int points_per_slice = 100; // How many points (on average) are in each slice of the 'pearl-necklace' sausage-length measurer
	double colloids[2][3]; // xyz positions of the two colloids
}

void invalid_colloid_info(void){
	error()<<"Invalid colloid information"<<std::endl;
	exit(EXIT_FAILURE);
}


void set_params(char *filename){
	std::ifstream paramfile (filename);
	if (!paramfile.is_open()){
		error()<<"Error opening file " << filename << std::endl;
		exit(EXIT_FAILURE);
	}

	std::stringstream buffer;
	buffer << paramfile.rdbuf();
	// Need silly c_str shenanigans to work with cJSON_Minify, which works at char* level
	char* c_str = new char[buffer.str().size()+1];
	std::strcpy(c_str,buffer.str().c_str());
	cJSON_Minify(c_str);
	cJSON *root=cJSON_Parse(c_str);
	delete [] c_str;


	if (!root){
		error()<<"Error reading param file before:" << std::endl << cJSON_GetErrorPtr() << std::endl;
		exit(EXIT_FAILURE);
	}

	/** @todo
	 * Could we do this nicely, with a dictionary or something? Types are a pain though.
	 * Could maybe use a union with operators? see http://www.dreamincode.net/forums/topic/129708-c-map-with-multiple-data-types/
	 * Would need to store the type along with default value.
	 * Or null pointers?
	 */

	/* Specials */
	// verbosity
	if (cJSON_GetObjectItem(root,"verbosity")){
			verbosity = (verbosityLevel) cJSON_GetObjectItem(root,"verbosity")->valueint;
			std::cout<<"verbosityLevel: ";
			switch(verbosity){
				case 0: std::cout<<"ERROR"; break;
				case 1: std::cout<<"WARNING"; break;
				case 2: std::cout<<"INFO"; break;
				case 3: std::cout<<"VERBOSE"; break;
				case 4: std::cout<<"DEBUG"; break;
			}
			std::cout<<std::endl;
	}

	// Colloids
	if (!cJSON_GetObjectItem(root,"colloids")){
		error()<<"No colloid information found in param file"<<std::endl;
		exit(EXIT_FAILURE);
	}

	cJSON *first = cJSON_GetObjectItem(root,"colloids")->child;
	if ( !first || !first->child || !first->child->next || !first->child->next->next )
		invalid_colloid_info();
	cJSON *second = first->next;
	if ( !second || !second->child || !second->child->next || !second->child->next->next )
		invalid_colloid_info();

	colloids[0][0]=first->child->valuedouble;
	colloids[0][1]=first->child->next->valuedouble;
	colloids[0][2]=first->child->next->next->valuedouble;
	colloids[1][0]=second->child->valuedouble;
	colloids[1][1]=second->child->next->valuedouble;
	colloids[1][2]=second->child->next->next->valuedouble;

	debug()<<"Colloid positions:"<<std::endl;
	for (int i=0; i<2; i++){
		for (int j=0; j<3; j++){
			debug()<<colloids[i][j]<<std::endl;
		}
	debug()<<std::endl;
	}



	/* Doubles */
	// threshold_high
	if (cJSON_GetObjectItem(root,"threshold_high")){
			threshold_high = cJSON_GetObjectItem(root,"threshold_high")->valuedouble;
			info()<<"threshold_high "<<threshold_high<<std::endl;
	}

	// threshold_low
	if (cJSON_GetObjectItem(root,"threshold_low")){
			threshold_low = cJSON_GetObjectItem(root,"threshold_low")->valuedouble;
			info()<<"threshold_low "<<threshold_low<<std::endl;
	}

	// silent_ignore_size
	if (cJSON_GetObjectItem(root,"silent_ignore_size")){
			silent_ignore_size = cJSON_GetObjectItem(root,"silent_ignore_size")->valuedouble;
			info()<<"silent_ignore_size "<<silent_ignore_size<<std::endl;
	}

	// min_sausage_size
	if (cJSON_GetObjectItem(root,"min_sausage_size")){
			min_sausage_size = cJSON_GetObjectItem(root,"min_sausage_size")->valuedouble;
			info()<<"min_sausage_size "<<min_sausage_size<<std::endl;
	}


	/* Integers */
	// points_per_slice
	if (cJSON_GetObjectItem(root,"points_per_slice")){
			points_per_slice = cJSON_GetObjectItem(root,"points_per_slice")->valuedouble;
			info()<<"points_per_slice "<<points_per_slice<<std::endl;
	}

}
