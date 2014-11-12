/** @file
 * This file holds all input/output functions.
 *
 * @todo
 * Eventually we should have .vtk and .gz capability,
 * and be able to distinguish dynamically which to use.
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "point.h"
#include "io.h"
#include "params.h"

using namespace std;

/** A null buffer and stream, analogous to /dev/null */
NullBuffer nullBuffer;
std::ostream nullStream(&nullBuffer);

/** Functions returning std::cout only at certain verbosity levels, else nullStream */
std::ostream &error(void) {
	if (params::verbosity >= ERROR){ return std::cout; } else{ return nullStream; }
}

std::ostream &warning(void) {
	if (params::verbosity >= WARNING){ return std::cout; } else{ return nullStream; }
}

std::ostream &brief(void) {
	if (params::verbosity >= BRIEF){ return std::cout; } else{ return nullStream; }
}

std::ostream &info(void) {
	if (params::verbosity >= INFO){ return std::cout; } else{ return nullStream; }
}

std::ostream &verbose(void) {
	if (params::verbosity >= VERBOSE){ return std::cout; } else{ return nullStream; }
}

std::ostream &debug(void) {
	if (params::verbosity >= DEBUG){ return std::cout; } else{ return nullStream; }
}



/**
 * Reads in data from inputFile with the format:
 * "x,y,z,cl,cp,cs", with one point per line.
 * Each point is appended to PointsVector, which is passed in by reference.
 * Each point's sausageID is returned undefined.
 */
void read_xyzclcpcs(std::ifstream& inputFile, vector<Point> &allPointsVector){
	string line;
	size_t iPoint=0;
	while ( getline (inputFile, line) ){
		Point p;
		sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf",
			&p.x,&p.y,&p.z,&p.cl,&p.cp,&p.cs);
		p.allPointsIndex=iPoint++;
		try{
			allPointsVector.push_back(p);
		} catch (std::bad_alloc &e){
			size_t pointsRead=allPointsVector.size();
			allPointsVector.clear();
			error()<<"Ran out of memory after reading in "<<pointsRead<<" points."<<std::endl;
			error()<<"Counting points left..."<<std::flush;
			size_t pointsLeft=0;
			while ( getline (inputFile, line) ){ pointsLeft++;}
			error()<<pointsLeft<<std::endl;
			error()<<"i.e. we could only allocate for "<<100*float(pointsRead)/float(pointsRead+pointsLeft)
				<<"% of the points in the file, increase memory available."<<std::endl;
			exit(EXIT_FAILURE);
		}

	}
}

