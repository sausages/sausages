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
#include "main.h"
#include "io.h"

#define INITPOINTARRAYSIZE 1000

using namespace std;


/**
 * Reads in data from inputFile with the format:
 * "x,y,z,cl,cp,cs", with one point per line.
 * Each point is appended to PointsArray, which is passed in by reference.
 * Each point's sausageID is returned undefined.
 */
void read_xyzclcpcs(std::ifstream& inputFile, Point** ptr_to_allPointsArray){
	Point *tmpArray=(Point*)calloc((size_t) INITPOINTARRAYSIZE,sizeof(Point));
	string line;
	Point p;
	size_t pointIndex = 0;
	while ( getline (inputFile, line) ){
		sscanf(line.c_str(),"%lf %lf %lf %lf %lf %lf",
			&p.x,&p.y,&p.z,&p.cl,&p.cp,&p.cs);
		cout << p.z << endl ;
		tmpArray[pointIndex]=p;
		pointIndex++;
	}
	*ptr_to_allPointsArray=tmpArray;
}
