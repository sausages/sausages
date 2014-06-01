#ifndef POINT_H
#define POINT_H

#include <vector>

struct Point
{
	double x;  ///< cartesian x-coord of point
	double y;  ///< cartesian y-coord of point
	double z;  ///< cartesian z-coord of point
	double cl; ///< raw cl value of point
	double cp; ///< raw cp value of point
	double cs; ///< raw cl value of point

	/** 0 for cl<threshold, 1 for unsorted.
	 * When sorted, each distinct sausage has a different number >=2.
	 * Each point is assigned to one of these sausages.
	 */
	int sausageID;

	/** Pointers to neighbours in 3D grid
	 * x:+right/-left ; y:+up/-down ; z:+forward/-back
	 * Don't need up-left etc. as can do up->left
	 */
	Point* left; 
	Point* right;
	Point* up;
	Point* down;
	Point* forward;
	Point* back;
};

void neighLink_xyzclcpcs(std::vector<Point> &allPoints);

void printAllNeighs(const std::vector<Point> &allPoints);

#endif // POINT_H
