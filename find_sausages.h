#ifndef FIND_SAUSAGES_H
#define FIND_SAUSAGES_H

#include <vector>
#include "point.h"

void threshold(std::vector<Point> &allPoints, double threshold_level);
std::vector<int> count_sausages(const std::vector<Point> &allPoints);
void flood_fill(std::vector<Point> &allPoints);

#endif // FIND_SAUSAGES_H
