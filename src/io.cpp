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
#include "main.h"
#include "maths.h"
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

/** Code relating to the 'Brief' file, standardised output info in a different file to cout */
bool briefIsInitialised = false;
std::ofstream &brief(void) { if (!briefIsInitialised) { initialiseBrief(); } return briefFile; }
std::ofstream briefFile;

void initialiseBrief(void){
	if (params::brief_filename.empty()){
		std::cerr << "Cannot use 'brief' file until after params have been initialised" << std::endl;
		exit(EXIT_FAILURE);
	}
	briefFile.open(params::brief_filename.c_str());
	briefIsInitialised=true;
	return;
}

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
		ifstream infile (inputFileName.c_str());
		if (!infile.is_open()){
			cerr<<"Error opening file "<<inputFileName<<endl;
			exit(EXIT_FAILURE);
		}
		istream& inStream = infile;

		if (inputFileName.rfind(".diot")==inputFileName.size()-5){
			if (params::colloidsInParamFile){
				cerr << "Cannot have colloid positions in param file when using .diot format" << endl;
				exit(EXIT_FAILURE);
			}
			return read_diot(inStream,allPointsVector);
		} else {
			if (!params::colloidsInParamFile){
				cerr << "Must have colloid positions in param file when not using .diot format" << endl;
				exit(EXIT_FAILURE);
			}
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
		/** TODO this is inefficient, should maybe only happen when being careful */
		int wordcount = 0;
	        stringstream ss( line );
		string word;
		while( ss >> word ) ++wordcount;
		if (wordcount != 6 && wordcount != 10){
			warning() << "Seem to be " << wordcount
				<< " words on line " << iPoint+1
				<< ", this is unexpected (expect 6 or 10)" <<endl;
		}

		Point p;
		sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf",
			&p.x,&p.y,&p.z,&p.cl,&p.cp,&p.cs);
		p.allPointsIndex=iPoint++;

		try{
			allPoints.push_back(p);
		} catch (std::bad_alloc &e){
			size_t pointsRead=allPoints.size();
			allPoints.clear();
			error()<<"Ran out of memory after reading in "<<pointsRead<<" points."<<std::endl;
			error()<<"Counting points left..."<<std::flush;
			size_t pointsLeft=0;
			while ( getline (input, line) ){ pointsLeft++;}
			error()<<pointsLeft<<std::endl;
			error()<<"i.e. we could only allocate for "<<100*float(pointsRead)/float(pointsRead+pointsLeft)
				<<"% of the points in the file, increase memory available."<<std::endl;
			exit(EXIT_FAILURE);
		}
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


/**
 * Reads in data from inputSteam with the .diot format.
 * Each point is appended to PointsVector, which is passed in by reference.
 */
int read_diot(std::istream &input, std::vector<Point> &allPoints){
	if (params::colloidsInParamFile){
		error()<<"Cannot have colloid info in the param file when using .diot format"<<std::endl;
		exit(EXIT_FAILURE);
	}
	string version;
	int numVoxels[3];
	float voxelSize[3];
	float lowBounds[3];
	string line,buffer,scanningMode;
	// Lines can appear in any order, so make sure we've seen them
	bool bVersion=0,bNumVoxels=0,bVoxelSize=0,bLowBounds=0,bBeginClCpCs=0;
	// Read in model information
	while ( getline (input, line) ){
		line = line.substr(0,line.find("#")); // Remove comments (anything after #)
		line.erase(0,line.find_first_not_of(" \t")); // Remove leading whitespace
		if (string::npos != line.find_last_not_of(" \t")) line.erase(line.find_last_not_of(" \t")+1); // Remove trailing whitespace
		if (line=="") continue; // Skip blank lines

		stringstream ss(line);

		if (line.find("version")!=string::npos){
			ss>>buffer;
			ss>>version;
			bVersion=1;
		} else if (line.find("numVoxels")!=string::npos){
			ss>>buffer;
			ss>>numVoxels[0];
			ss>>numVoxels[1];
			ss>>numVoxels[2];
			bNumVoxels=1;
		} else if (line.find("voxelSize")!=string::npos){
			ss>>buffer;
			ss>>voxelSize[0];
			ss>>voxelSize[1];
			ss>>voxelSize[2];
			bVoxelSize=1;
		} else if (line.find("lowBounds")!=string::npos){
			ss>>buffer;
			ss>>lowBounds[0];
			ss>>lowBounds[1];
			ss>>lowBounds[2];
			bLowBounds=1;
		} else if (line.find("colloidPos")!=string::npos){
			ss>>buffer;
			size_t colloidNum;
			ss>>colloidNum;
			if (colloidNum!=model::colloidPos.size()+1){
				error()<<"Colloids must be declared in order, aborting"<<endl;
				exit(EXIT_FAILURE);
			}
			vector3d newColloid;
			ss>>newColloid.x;
			ss>>newColloid.y;
			ss>>newColloid.z;
			model::colloidPos.push_back(newColloid);
		} else if (line.find("colloidRad")!=string::npos){
			ss>>buffer;
			size_t colloidNum;
			ss>>colloidNum;
			if (colloidNum>model::colloidPos.size()){
				error()<<"Colloids properties must be declared after its position, aborting"<<endl;
				exit(EXIT_FAILURE);
			}
			double colloidRad;
			ss>>colloidRad;
			model::colloidRadii.push_back(colloidRad);
		} else if (line.find("beginClCpCs")!=string::npos){
			bBeginClCpCs=1;
			ss>>buffer;
			ss>>scanningMode;
			break;
		}
	}
	// Did we read everything we needed to?
	if (!(bVersion && bNumVoxels && bVoxelSize && bLowBounds && bBeginClCpCs)) {
		error()<<"Invalid .diot format"<<endl;
		exit(EXIT_FAILURE);
	}
	if (model::colloidPos.size()!=2){
		error()<<"Can currently only handle systems with 2 colloids, sorry"<<endl;
		exit(EXIT_FAILURE);
	}

	// Read in points, only saving threshold-passing ones, linking them properly
	int points_read=0; // !=allPoints.size();
	int ix=0,iy=0,iz=0; // Point coordinate in point-space
	double cl,cp,cs;

	string zyxInc ("zyxInc");
	if (scanningMode!=zyxInc){ // Lovely, lovely C++ strings <3
		error()<<"Currently can only handle zyxInc scanning mode in diot files."<<endl;
		exit(EXIT_FAILURE);
	}

	int iPoint=0;
	vector<int> pointMap; // Has a sensible structure, each element point to an allPoints element (or -1 if cl<threshold)
	while ( getline (input, line) ){
		points_read++;
		// Assuming BeginClCpCs is in zyxInc
		iz++;
		if (iz==numVoxels[2]){
			iz=0; iy++;
		}
		if (iy==numVoxels[1]){
			iy=0; ix++;
		}

		stringstream ss(line);
		ss>>cl;
		ss>>cp;
		ss>>cs;

		// Only save points that pass threshold condition
		if (cl < params::threshold){
			Point p;
			p.cl=cl; p.cs=cs; p.cp=cp;
			p.x = ix*voxelSize[0] + lowBounds[0];
			p.y = iy*voxelSize[1] + lowBounds[1];
			p.z = iz*voxelSize[2] + lowBounds[2];
			p.allPointsIndex = iPoint++;
			p.isInASausage=true;
			allPoints.push_back(p);
			pointMap.push_back(p.allPointsIndex);
		} else {
			pointMap.push_back(-1);
		}
	}

	/* Links each point in allPoints to its 6 nearest-neighbours, using pointMap.
	 * pointMap has an index-structure we can use, and each element point to an allPoints element
	 * x:+right/-left ; y:+up/-down ; z:+forward/-back
	 * Assumes scanningMode = zyxInc
	 */

	for (size_t i=0; i<pointMap.size(); i++){
		if (pointMap[i]<0) continue; // These are points w. cl<threshold

		int iCurr=pointMap[i]; // allPoint index
		int Nx = numVoxels[0], Ny = numVoxels[1], Nz = numVoxels[2];

		// Pointer to self, pretty sure this is irrelevant, but
		// I currently need it for the flood-fill (iterators are copies)
		allPoints[iCurr].self = &(allPoints[iCurr]);

		// If not right-most
		if (((i/(Nz*Ny))+1)%Nx != 0){
			if (pointMap[i+Nz*Ny] != -1){ // If right-of-me passed thresholding
				int iRight = pointMap[i+Nz*Ny];
				allPoints[iCurr].right  = &(allPoints[iRight]);
				allPoints[iRight].left = &(allPoints[iCurr]);
			}
		}
		// If not upper-most
		if (((i/Nz)+1)%Ny != 0){
			if (pointMap[i+Nz] != -1){ // If above-me passed thresholding
				int iUp = pointMap[i+Nz];
				allPoints[iCurr].up   = &(allPoints[iUp]);
				allPoints[iUp].down = &(allPoints[iCurr]);
			}
		}
		// If not forward-most
		if ((i+1)%Nz != 0){
			if (pointMap[i+1] != -1){ // If forward-of-me passed thresholding
				int iForward = pointMap[i+1];
				allPoints[iCurr].forward = &(allPoints[iForward]);
				allPoints[iForward].back   = &(allPoints[iCurr]);
			}
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
	return points_read;
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
	int points_read;
	if (filename.rfind(".diot")==filename.size()-5){
		points_read = read_diot(iss,allPointsVector);
	} else {
		points_read = read_xyzclcpcs(iss,allPointsVector);
	}

	mz_free(p);
	mz_zip_reader_end(&zip_archive);

	return points_read;

}
