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
#include <cmath>
#include <vector>
#include "point.h"
#include "sausages.h"
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
int read_input(std::string inputFileName, std::vector<Point> &allPointsVector){
	// If inputFileName ends with ".zip"
	if (inputFileName.rfind(".zip")==inputFileName.size()-4){
		debug()<<"Looks like a zip file"<<endl;
		return read_zipped(inputFileName,allPointsVector);
	} else {
		debug()<<"Not a zip file"<<endl;
		ifstream infile (inputFileName);
		if (!infile.is_open()){
			cerr<<"Error opening file "<<inputFileName<<endl;
			exit(EXIT_FAILURE);
		}
		istream& inStream = infile;

		if (inputFileName.rfind(".diot")==inputFileName.size()-5){
			return 0;
	//		return read_diot(inStream,allPointsVector);
		} else {
			return read_xyzclcpcs(inStream,allPointsVector);
		}
	}

}


/**
 * Reads in data from inputSteam with the format:
 * "x,y,z,cl,cp,cs", with one point per line.
 * Each point is appended to PointsVector, which is passed in by reference.
 */
int read_xyzclcpcs(std::istream &input, std::vector<Point> &allPoints){
	string line;
	size_t iPoint=0;
	while ( getline (input, line) ){
		Point p;
		sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf",
			&p.x,&p.y,&p.z,&p.cl,&p.cp,&p.cs);
		p.allPointsIndex=iPoint++;
		allPoints.push_back(p);
	}

	/* Links each point in allPoints to its 6 nearest-neighbours.
	 * x:+right/-left ; y:+up/-down ; z:+forward/-back
	 * @warning Assumes that points are listed in a specific order:
	 * @warning   z varies first (increasing), then y, then x.
	 * @warning Also assumes that grid is cubic, i.e. N=len(x)=len(y)=len(z)
	 */
	unsigned int N = round(pow(allPoints.size(),1.0/3.0));///< side-length of cube
	// Don't use iterators, indices are important here
	for (size_t i=0; i<allPoints.size(); i++){
		// Pointer to self, pretty sure this is irrelevant, but
		// I currently need it for the flood-fill (iterators are copies)
		allPoints[i].self = &(allPoints[i]);
		// If not right-most
		if (((i/(N*N))+1)%N != 0){
			allPoints[i].right  = &(allPoints[i+N*N]);
			allPoints[i+N*N].left = &(allPoints[i]);
		}
		// If not upper-most
		if (((i/N)+1)%N != 0){
			allPoints[i].up   = &(allPoints[i+N]);
			allPoints[i+N].down = &(allPoints[i]);
		}
		// If not forward-most
		if ((i+1)%N != 0){
			allPoints[i].forward = &(allPoints[i+1]);
			allPoints[i+1].back   = &(allPoints[i]);
		}
	}

	// Update neighbours array
	for (size_t i=0; i<allPoints.size(); i++){
		allPoints[i].neighbours[0]=allPoints[i].left;
		allPoints[i].neighbours[1]=allPoints[i].right;
		allPoints[i].neighbours[2]=allPoints[i].up;
		allPoints[i].neighbours[3]=allPoints[i].down;
		allPoints[i].neighbours[4]=allPoints[i].forward;
		allPoints[i].neighbours[5]=allPoints[i].back;
	}

	// Check Pixel size
	if (abs(abs(allPoints[0].z-allPoints[1].z)-params::pixel_size)>params::epsilon){
		warning()<<"WARNING: Pixel size appears to be "<<abs(allPoints[0].z-allPoints[1].z)<<
			" but params file (or default) is set to "<<params::pixel_size<<"."<<endl<<
			"Overwriting with new pixel size."<<endl;
		params::pixel_size=abs(abs(allPoints[0].z-allPoints[1].z)-params::pixel_size);
	}

	info()<<"Read in "<<allPoints.size()<<" points."<<endl;
	info()<<"Thresholding, cl<"<<params::threshold<<" belong to a sausage"<<endl;
	int numInASausage = threshold(allPoints);
	info()<<numInASausage<<" points left after thresholding."<<endl;
	return numInASausage;
}


int read_zipped(std::string inputArchiveFileName, std::vector<Point> &allPointsVector){
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
	string filename(file_stat.m_filename);
	if (filename.rfind(".diot")==filename.size()-5){
		return 0;
		//return read_diot(iss,allPointsVector);
	} else {
		return read_xyzclcpcs(iss,allPointsVector);
	}

	mz_free(p);
	mz_zip_reader_end(&zip_archive);

}
