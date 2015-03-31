#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cstring>
#include "cJSON/cJSON.h"
#include "io.h"
#include "main.h"
#include "maths.h"
#include "params.h"

using namespace params;

// Set defaults
namespace params{
	verbosityLevel verbosity = INFO;
	std::string brief_filename = ""; // 'brief' file is for standardised output, in a different file to cout
	int brief_version = 1 ; // each brief version is standardised
	std::string thresholded_filename = ""; // 'thresholded' file is for output of points below the threshold, in a different file to cout
	std::string sausage_filename = ""; // 'sausage' file is for output of points in sausages, one file per sausage
	double threshold = 0.04; // c_l threshold for includion into a sausage. All different defect structures should be distinguishable (i.e. no ambigious blobs)
	double silent_ignore_size = 0.01; // If a sausage is smaller than this fraction of all points below the threshold, silently ignore it
	double min_sausage_size = 0.1; // A sausage is only 'significant' if it is larger than this fraction of all points below the threshold
	double pixel_size = 0.5; // Distance between nearest-neighbour points
	bool colloidsInParamFile = false; // Colloids shouldn't be both in param and input file
	double ratio_two_rings = 0.1; // Threshold ratio for which two relevant sausages are identified as two_rings
	double ratio_2nd_loop = 0.5; // Threshold ratio for which two relevant sausages are identified as 2nd_loop
	double flood_fill_classify_slice_size=4; // How many pixels wide should the regions in flood_fill_classify be?
	double epsilon=1.0e-10; // Some small number for comparison to zero
	double max_sausage_gap_size=15; // Maximum distance (in units of pixel-length) between endpoints that will be joined by join_endpoints()
	double min_sausage_gap_next_neigh_distance=20; // If there are multiple neighbouring endpoints within this radius (in units of pixels), we throw an error as it's too close to call between candidates
    double R_min_sphere = 0.1; // Minimum radius for sphere tracking algorithm
    double dR_sphere = 0.1; // Increament sphere by dR value in sphere tracking algorithm
    double R_gap = 1.0; // Small radius used in sphere tracking to determine our next point, distance between adjacent points ~ R_gap
}

void invalid_colloid_info(void){
	error()<<"Invalid colloid information"<<std::endl;
	exit(EXIT_FAILURE);
}


void set_params(char *filename, std::vector<vector3d>  &colloidPos){
	if (filename == NULL){
		// Default initialisation
		brief_filename=std::string("default.brief");
		return;
	}

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
				case 2: std::cout<<"BRIEF";break;
				case 3: std::cout<<"INFO"; break;
				case 4: std::cout<<"VERBOSE"; break;
				case 5: std::cout<<"DEBUG"; break;
			}
			std::cout<<std::endl;
	}

	// Brief filename
	if (cJSON_GetObjectItem(root,"brief_filename")){
		brief_filename=std::string(cJSON_GetObjectItem(root,"brief_filename")->valuestring);
	} else {
		brief_filename=std::string("default.brief");
	}
	info()<<"brief_filename "<<brief_filename<<std::endl;

	// points_filename
	if (cJSON_GetObjectItem(root,"thresholded_filename")){
			thresholded_filename = cJSON_GetObjectItem(root,"thresholded_filename")->valuestring;
			info()<<"thresholded_filename "<<thresholded_filename<<std::endl;
	}

	// sausage_filename
	if (cJSON_GetObjectItem(root,"sausage_filename")){
			sausage_filename = cJSON_GetObjectItem(root,"sausage_filename")->valuestring;
			info()<<"sausage_filename "<<sausage_filename<<std::endl;
	}


	// Colloids, we can only deal with two in a param file (this is deprecated, pre-DIOT)
	if (!cJSON_GetObjectItem(root,"colloids")){
		colloidsInParamFile=false;
	} else {
		colloidsInParamFile=true;
		cJSON *first = cJSON_GetObjectItem(root,"colloids")->child;
		if ( !first || !first->child || !first->child->next || !first->child->next->next )
			invalid_colloid_info();
		cJSON *second = first->next;
		if ( !second || !second->child || !second->child->next || !second->child->next->next )
			invalid_colloid_info();

		for (int i=0; i<2; i++){
			cJSON *currColloid = (i==0) ? first : second ;
			vector3d newColloid;
			newColloid.x=currColloid->child->valuedouble;
			newColloid.y=currColloid->child->next->valuedouble;
			newColloid.z=currColloid->child->next->next->valuedouble;
			colloidPos.push_back(newColloid);
		}

		debug()<<"Colloid positions:"<<std::endl;
		for (int i=0; i<2; i++){
			debug()<<colloidPos[i]<<std::endl;
		}
	}



	/* Doubles */
	// threshold
	if (cJSON_GetObjectItem(root,"threshold")){
			params::threshold= cJSON_GetObjectItem(root,"threshold")->valuedouble;	// Need to specify threshold to avoid conflict
			info()<<"threshold "<<params::threshold<<std::endl;			// with function in sausages.h, via main.h
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

	// pixel_size
	if (cJSON_GetObjectItem(root,"pixel_size")){
			pixel_size = cJSON_GetObjectItem(root,"pixel_size")->valuedouble;
			info()<<"pixel_size "<<pixel_size<<std::endl;
	}

	// flood_fill_classify_slice_size
	if (cJSON_GetObjectItem(root,"flood_fill_classify_slice_size")){
			flood_fill_classify_slice_size = cJSON_GetObjectItem(root,"flood_fill_classify_slice_size")->valuedouble;
			info()<<"flood_fill_classify_slice_size "<<flood_fill_classify_slice_size<<std::endl;
	}

	// flood_fill_classify_slice_size
	if (cJSON_GetObjectItem(root,"epsilon")){
			epsilon = cJSON_GetObjectItem(root,"epsilon")->valuedouble;
			info()<<"epsilon "<<epsilon<<std::endl;
	}


	/* Integers */
	if (cJSON_GetObjectItem(root,"brief_version")){
			brief_version = cJSON_GetObjectItem(root,"brief_version")->valuedouble;
			info()<<"brief_version "<<brief_version<<std::endl;
	}

	// ratio_two_rings
	if (cJSON_GetObjectItem(root,"ratio_two_rings")){
			ratio_two_rings = cJSON_GetObjectItem(root,"ratio_two_rings")->valuedouble;
			info()<<"ratio_two_rings "<<ratio_two_rings<<std::endl;
	}

	// ratio_2nd_loop
	if (cJSON_GetObjectItem(root,"ratio_2nd_loop")){
			ratio_2nd_loop = cJSON_GetObjectItem(root,"ratio_2nd_loop")->valuedouble;
			info()<<"ratio_2nd_loop "<<ratio_2nd_loop<<std::endl;
	}


	// Clean up
	cJSON_Delete(root);
	paramfile.close();

}
