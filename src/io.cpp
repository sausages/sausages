/** @file
 * This file holds all input/output functions.
 *
 * @todo
 * Eventually we should have .vtk and .gz capability,
 * and be able to distinguish dynamically which to use.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "point.h"
#include "io.h"
#include "params.h"
#include "miniz/miniz.c"

using namespace std;

/** A null buffer and stream, analogous to /dev/null */
NullBuffer nullBuffer;
std::ostream nullStream(&nullBuffer);

/** Functions returning std::cout only at certain verbosity levels, else nullStream */
std::ostream &error(void) { if (params::verbosity >= ERROR){ return std::cout; } else{ return nullStream; } }
std::ostream &warning(void) { if (params::verbosity >= WARNING){ return std::cout; } else{ return nullStream; } }
std::ostream &info(void) { if (params::verbosity >= INFO){ return std::cout; } else{ return nullStream; } }
std::ostream &verbose(void) { if (params::verbosity >= VERBOSE){ return std::cout; } else{ return nullStream; } }
std::ostream &debug(void) { if (params::verbosity >= DEBUG){ return std::cout; } else{ return nullStream; } }

/** Decide, based on filename extension, what file type the input is,
 *  and call the relevant function
 */
void read_file(std::string inputFileName, std::vector<Point> &allPointsVector){
	// If inputFileName ends with ".zip"
	if (inputFileName.rfind(".zip")==inputFileName.size()-4){
		debug()<<"Looks like a zip file"<<endl;
		read_zipped_xyzclcpcs(inputFileName,allPointsVector);
	} else {
		debug()<<"Not a zip file"<<endl;
		ifstream infile (inputFileName);
		if (!infile.is_open()){
			cerr<<"Error opening file "<<inputFileName<<endl;
			exit(EXIT_FAILURE);
		}
		istream& inStream = infile;
		read_xyzclcpcs(inStream,allPointsVector);
	}

}


/**
 * Reads in data from inputSteam with the format:
 * "x,y,z,cl,cp,cs", with one point per line.
 * Each point is appended to PointsVector, which is passed in by reference.
 */
void read_xyzclcpcs(std::istream &input, std::vector<Point> &allPointsVector){
	string line;
	size_t iPoint=0;
	while ( getline (input, line) ){
		Point p;
		sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf",
			&p.x,&p.y,&p.z,&p.cl,&p.cp,&p.cs);
		p.allPointsIndex=iPoint++;
		allPointsVector.push_back(p);
	}
}


void read_zipped_xyzclcpcs(std::string inputArchiveFileName, std::vector<Point> &allPointsVector){
	// Initialise
	mz_zip_archive zip_archive;
	memset(&zip_archive,0,sizeof(zip_archive));
	mz_bool status = mz_zip_reader_init_file(&zip_archive,inputArchiveFileName.c_str(),0);
	if (!status){
		error()<<"Couldn't initialise zip-reader on file "<<inputArchiveFileName<<endl;
		exit(EXIT_FAILURE);
	}

	// Check number of files inside archive, if there's more than one we don't know which to use, so quit
	if (mz_zip_reader_get_num_files(&zip_archive) != 1){
		error()<<"File "<<inputArchiveFileName<<" appears to contain more than a single file, aborting"<<endl;
		exit(EXIT_FAILURE);
	}

	// Get compressed file's details
	mz_zip_archive_file_stat file_stat;
	mz_zip_reader_file_stat(&zip_archive,0,&file_stat);
	debug()<<"Found file "<<file_stat.m_filename<<" in archive"<<endl;

	// Extract file
	size_t uncomp_size;
	void *p;
	p=mz_zip_reader_extract_file_to_heap(&zip_archive,file_stat.m_filename, &uncomp_size,0);
	if (!p){
		error()<<"Failed to extract file "<<file_stat.m_filename<<" from archive "<<inputArchiveFileName<<endl;
		exit(EXIT_FAILURE);
	}

	debug()<<"Successfully extracted file"<<endl;

	istringstream iss;
	iss.str((char*)p);
	read_xyzclcpcs(iss,allPointsVector);

	mz_free(p);
	mz_zip_reader_end(&zip_archive);

}
