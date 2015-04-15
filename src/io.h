#ifndef IO_H
#define IO_H

#include <iostream>
#include <fstream>
#include <vector>
#include "point.h"
#include "main.h"

void read_input(std::string inputFileName, Model &model);
void read_xyzclcpcs(std::istream& input, Model &model);
void read_zipped(std::string inputArchiveFileName, Model &model);
void read_diot(std::istream &input, Model &model);

void write_points(std::string filename, std::vector<Point> points);
void write_points(std::string filename, std::vector<Point*> points);
// Same as above, but print only the elements of 'points' whose index is in 'indices'
void write_points(std::string filename, std::vector<Point*> points, std::vector<int> indices);

class NullBuffer : public std::streambuf
{
	public:
		  int overflow(int c) { return c; }
};

extern bool briefIsInitialised;
extern std::ostream &brief(std::vector<int> versions);
extern std::ofstream briefFile;
void initialiseBrief(void);

/** Verious dummy streams, which will output to cout only when verbosityLevel is high enough */
extern std::ostream &error(void) ;
extern std::ostream &warning(void) ;
extern std::ostream &info(void) ;
extern std::ostream &verbose(void) ;
extern std::ostream &debug(void) ;

#endif // IO_H
