#ifndef IO_H
#define IO_H

#include <iostream>
#include <fstream>
#include <vector>
#include "point.h"

int read_input(std::string inputFileName, std::vector<Point> &allPointsVector);
int read_xyzclcpcs(std::istream& input, std::vector<Point> &allPointsVector);
int read_zipped(std::string inputArchiveFileName, std::vector<Point> &allPointsVector);
int read_diot(std::istream &input, std::vector<Point> &allPoints);

class NullBuffer : public std::streambuf
{
	public:
		  int overflow(int c) { return c; }
};

extern bool briefIsInitialised;
extern std::ofstream &brief(void);
extern std::ofstream briefFile;
void initialiseBrief(void);

/** Verious dummy streams, which will output to cout only when verbosityLevel is high enough */
extern std::ostream &error(void) ;
extern std::ostream &warning(void) ;
extern std::ostream &info(void) ;
extern std::ostream &verbose(void) ;
extern std::ostream &debug(void) ;

#endif // IO_H
