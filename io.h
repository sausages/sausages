#ifndef IO_H
#define IO_H

#include <iostream>
#include <fstream>
#include "main.h"

void read_xyzclcpcs(std::ifstream& inputFile, Point** ptr_to_allPointsArray);

#endif // IO_H
