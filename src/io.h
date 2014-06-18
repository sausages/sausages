#ifndef IO_H
#define IO_H

#include <iostream>
#include <fstream>
#include <vector>
#include "point.h"

void read_xyzclcpcs(std::ifstream& inputFile, std::vector<Point> &allPointsVector);

class NullBuffer : public std::streambuf
{
	public:
		  int overflow(int c) { return c; }
};

/** Verious dummy streams, which will output to cout only when verbosityLevel is high enough */
extern std::ostream &error(void) ;
extern std::ostream &warning(void) ;
extern std::ostream &info(void) ;
extern std::ostream &verbose(void) ;
extern std::ostream &debug(void) ;

#endif // IO_H
