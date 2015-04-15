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

/** Code relating to the 'Brief' file, standardised output info in a different file to cout
 * Usage:
 * brief({1,2,4}) << "This is printed in versions 1,2,4 only" << endl;
 */

bool briefIsInitialised = false;
std::ofstream briefFile;

std::ostream &brief(std::vector<int> versions) {
	if (!briefIsInitialised) initialiseBrief();
	if (elementOf(versions, params::brief_version)){
		return briefFile;
	} else {
		return nullStream;
	}
}

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
void read_input(std::string inputFileName, Model &model){
	// If inputFileName ends with ".zip"
	if (inputFileName.rfind(".zip")==inputFileName.size()-4){
		debug()<<"Looks like a zip file"<<endl;
		return read_zipped(inputFileName, model);
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
			read_diot(inStream, model);
		} else {
			if (!params::colloidsInParamFile){
				cerr << "Must have colloid positions in param file when not using .diot format" << endl;
				exit(EXIT_FAILURE);
			}
			read_xyzclcpcs(inStream, model);
		}
	}

	return;
}


/**
 * Reads in data from inputSteam with the format:
 * "x,y,z,cl,cp,cs", with one point per line.
 * Each point is appended to PointsVector, which is passed in by reference.
 */
void read_xyzclcpcs(std::istream &input, Model &model){
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

		Point &p = *(new Point ());
		sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf",
			&p.x,&p.y,&p.z,&p.cl,&p.cp,&p.cs);
		p.allPointsIndex=iPoint++;

		try{
			model.allPoints.push_back(&p);
		} catch (std::bad_alloc &e){
			size_t pointsRead=model.allPoints.size();
			for (Point* pp : model.allPoints) {delete pp;}
			model.allPoints.clear();
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
	model.total_points = model.allPoints.size();

	/* Links each point in allPoints to its 6 nearest-neighbours.
	 * x:+right/-left ; y:+up/-down ; z:+forward/-back
	 * @warning Assumes that points are listed in a specific order:
	 * @warning   z varies first (increasing), then y, then x.
	 * @warning Also assumes that grid is cubic, i.e. N=len(x)=len(y)=len(z)
	 */
	unsigned int N = round(pow(model.allPoints.size(),1.0/3.0));///< side-length of cube
	// Don't use iterators, indices are important here
	for (size_t i=0; i<model.allPoints.size(); i++){
		// If not right-most
		if (((i/(N*N))+1)%N != 0){
			model.allPoints[i]->right  = model.allPoints[i+N*N];
			model.allPoints[i+N*N]->left = model.allPoints[i];
		}
		// If not upper-most
		if (((i/N)+1)%N != 0){
			model.allPoints[i]->up   = model.allPoints[i+N];
			model.allPoints[i+N]->down = model.allPoints[i];
		}
		// If not forward-most
		if ((i+1)%N != 0){
			model.allPoints[i]->forward = model.allPoints[i+1];
			model.allPoints[i+1]->back   = model.allPoints[i];
		}
	}

	// Update neighbours array
	for (size_t i=0; i<model.allPoints.size(); i++){
		model.allPoints[i]->neighbours[0]=model.allPoints[i]->left;
		model.allPoints[i]->neighbours[1]=model.allPoints[i]->right;
		model.allPoints[i]->neighbours[2]=model.allPoints[i]->up;
		model.allPoints[i]->neighbours[3]=model.allPoints[i]->down;
		model.allPoints[i]->neighbours[4]=model.allPoints[i]->forward;
		model.allPoints[i]->neighbours[5]=model.allPoints[i]->back;
	}

	// Check Pixel size
	if (abs(abs(model.allPoints[0]->z-model.allPoints[1]->z)-params::pixel_size)>params::epsilon){
		warning()<<"WARNING: Pixel size appears to be "<<abs(model.allPoints[0]->z - model.allPoints[1]->z)<<
			" but params file (or default) is set to "<<params::pixel_size<<"."<<endl<<
			"WARNING: Overwriting with new pixel size."<<endl;
		params::pixel_size=abs(model.allPoints[0]->z - model.allPoints[1]->z);
	}

	// Points are in a sausge if their cl value is < threshold
	debug() << "begin threshold" << endl << flush;

	for (Point* p : model.allPoints){
		if (p->cl < params::threshold){
			p->isInASausage=true;
			model.threshold_points.push_back(p->allPointsIndex);
		}
	}
	debug() << "end threshold" << endl << flush;

	return;
}


/**
 * Reads in data from inputSteam with the .diot format.
 * Each point is appended to PointsVector, which is passed in by reference.
 */
void read_diot(std::istream &input, Model &model){
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
			if (colloidNum!=model.colloidPos.size()+1){
				error()<<"Colloids must be declared in order, aborting"<<endl;
				exit(EXIT_FAILURE);
			}
			Vector3d newColloid (0,0,0);
			ss>>newColloid.x;
			ss>>newColloid.y;
			ss>>newColloid.z;
			model.colloidPos.push_back(newColloid);
		} else if (line.find("colloidRad")!=string::npos){
			ss>>buffer;
			size_t colloidNum;
			ss>>colloidNum;
			if (colloidNum>model.colloidPos.size()){
				error()<<"Colloids properties must be declared after its position, aborting"<<endl;
				exit(EXIT_FAILURE);
			}
			double colloidRad;
			ss>>colloidRad;
			model.colloidRadii.push_back(colloidRad);
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
	if (model.colloidPos.size()!=2){
		error()<<"Can currently only handle systems with 2 colloids, sorry"<<endl;
		exit(EXIT_FAILURE);
	}

	// Read in points, only saving threshold-passing ones, linking them properly
	model.total_points=0; // !=allPoints.size();
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
		model.total_points++;
		// Assuming BeginClCpCs is in zyxInc
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
			Point &p = *(new Point () );
			p.cl=cl; p.cs=cs; p.cp=cp;
			p.x = ix*voxelSize[0] + lowBounds[0];
			p.y = iy*voxelSize[1] + lowBounds[1];
			p.z = iz*voxelSize[2] + lowBounds[2];
			p.allPointsIndex = iPoint++;
			p.isInASausage=true;
			model.allPoints.push_back(&p);
			model.threshold_points.push_back(p.allPointsIndex);
			pointMap.push_back(p.allPointsIndex);
		} else {
			pointMap.push_back(-1);
		}

		// Assuming BeginClCpCs is in zyxInc
		iz++;
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

		// If not right-most
		if (((i/(Nz*Ny))+1)%Nx != 0){
			if (pointMap[i+Nz*Ny] != -1){ // If right-of-me passed thresholding
				int iRight = pointMap[i+Nz*Ny];
				model.allPoints[iCurr]->right  = model.allPoints[iRight];
				model.allPoints[iRight]->left = model.allPoints[iCurr];
			}
		}
		// If not upper-most
		if (((i/Nz)+1)%Ny != 0){
			if (pointMap[i+Nz] != -1){ // If above-me passed thresholding
				int iUp = pointMap[i+Nz];
				model.allPoints[iCurr]->up   = model.allPoints[iUp];
				model.allPoints[iUp]->down = model.allPoints[iCurr];
			}
		}
		// If not forward-most
		if ((i+1)%Nz != 0){
			if (pointMap[i+1] != -1){ // If forward-of-me passed thresholding
				int iForward = pointMap[i+1];
				model.allPoints[iCurr]->forward = model.allPoints[iForward];
				model.allPoints[iForward]->back   = model.allPoints[iCurr];
			}
		}
	}

	// Update neighbours array
	for (size_t i=0; i<model.allPoints.size(); i++){
		model.allPoints[i]->neighbours[0]=model.allPoints[i]->left;
		model.allPoints[i]->neighbours[1]=model.allPoints[i]->right;
		model.allPoints[i]->neighbours[2]=model.allPoints[i]->up;
		model.allPoints[i]->neighbours[3]=model.allPoints[i]->down;
		model.allPoints[i]->neighbours[4]=model.allPoints[i]->forward;
		model.allPoints[i]->neighbours[5]=model.allPoints[i]->back;
	}

	return;
}


void read_zipped(std::string inputArchiveFileName, Model &model){
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
		read_diot(iss, model);
	} else {
		read_xyzclcpcs(iss, model);
	}

	mz_free(p);
	mz_zip_reader_end(&zip_archive);

	return;

}


/** Gnuplot (splot) is the nicest way to view lists of vectors,
 *  and it prefers 'x y z' per line to Point-stream-overload of '(x,y,z)'
 */
void write_points(std::string filename, std::vector<Point> points){
	debug() << "Printing to file " << filename << endl;
	std::ofstream outfile (filename);
	for (Point iter: points){
		outfile << iter.x << ' ' << iter.y << ' ' << iter.z << endl;
	}
	outfile.close();
}
// Same as above, but with vectors-of-pointers
void write_points(std::string filename, std::vector<Point*> points){
	debug() << "Printing to file " << filename << endl;
	std::ofstream outfile (filename);
	for (Point* iter: points){
		outfile << iter->x << ' ' << iter->y << ' ' << iter->z << endl;
	}
	outfile.close();
}
// Same as above, but print only the elements of 'points' whose index is in 'indices'
void write_points(std::string filename, std::vector<Point*> points, std::vector<int> indices){
	debug() << "Printing to file " << filename << endl;
	std::ofstream outfile (filename);
	for (int i : indices){
		Point* iter = points[i];
		outfile << iter->x << ' ' << iter->y << ' ' << iter->z << endl;
	}
	outfile.close();
}
