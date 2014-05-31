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
#include "main.h"
#include "io.h"

using namespace std;


/**
 * Reads in data from inputFile with the format:
 * "x,y,z,cl,cp,cs", with one point per line.
 * Each point is appended to PointsVector, which is passed in by reference.
 * Each point's sausageID is returned undefined.
 */
void read_xyzclcpcs(std::ifstream& inputFile, vector<Point> &allPointsVector){
	string line;
	Point p;
	while ( getline (inputFile, line) ){
		sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf",
			&p.x,&p.y,&p.z,&p.cl,&p.cp,&p.cs);
		allPointsVector.push_back(p);
	}
}
