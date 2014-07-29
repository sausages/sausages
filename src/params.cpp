#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "cJSON/cJSON.h"
#include "io.h"
#include "params.h"

using namespace params;

// Set defaults (see params.h for descriptions)
namespace params{
	const verbosityLevel verbosity = INFO;
	const double threshold_level = 0.12;
	const double silent_ignore_size = 0.01;
	const double min_sausage_size = 0.1;
	const int pointsPerSlice = 100;
}

void set_params(char *filename){
	std::ifstream paramfile (filename);
	if (!paramfile.is_open()){
		error()<<"Error opening file " << filename << std::endl;
		exit(EXIT_FAILURE);
	}

	std::stringstream buffer;
	buffer << paramfile.rdbuf();
	cJSON *json=cJSON_Parse(buffer.str().c_str());

	if (!json){
		error()<<"Error before: [" << cJSON_GetErrorPtr() << "]" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::cout << cJSON_GetObjectItem(json,"myint")->valueint << std::endl;
}
