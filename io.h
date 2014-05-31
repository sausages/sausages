#ifndef IO_H
#define IO_H

#include <iostream>
#include <fstream>
#include <vector>
#include "main.h"

void read_xyzclcpcs(std::ifstream& inputFile, std::vector<Point> &allPointsVector);

#endif // IO_H
