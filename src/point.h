#ifndef POINT_H
#define POINT_H

#include <vector>
#include "maths.h"

class Point
{
	public:

	double x;  ///< cartesian x-coord of point
	double y;  ///< cartesian y-coord of point
	double z;  ///< cartesian z-coord of point
	double cl; ///< raw cl value of point
	double cp; ///< raw cp value of point
	double cs; ///< raw cl value of point

	operator Vector3d(); ///< Cast to Vector3d returns only x,y,z coords

	int sausageID; ///< Which sausage the point is in (-1 is unsorted)
	bool isInASausage;

	// The index of this point in the allPoints array (so we don't have to search)
	size_t allPointsIndex;

	// The index of this point in its sauasages points array (so we don't have to search)
	size_t sausagePointsIndex;

	/** Pointers to neighbours in 3D grid
	 * x:+right/-left ; y:+up/-down ; z:+forward/-back
	 * Don't need up-left etc. as can do up->left
	 * These should be referenced only during setup, in general
	 * loop over the neighbours array.
	 */
	Point* left; 
	Point* right;
	Point* up;
	Point* down;
	Point* forward;
	Point* back;

	/** Array of {l,r,u,d,f,b}
	 * Advantage of array is possibility of looping over it
	 */
	Point* neighbours[6];

	// Assume if our allPointsIndex is the same then we are the same
	// Currently needed for the find algorithm in flood-fill
	bool operator == (const Point &Ref) const {
		return(this->allPointsIndex == Ref.allPointsIndex);
	}
	// Sort points (for uniqueness) by value of allPointsIndex
	bool operator < (const Point &Ref) const {
		return(this->allPointsIndex < Ref.allPointsIndex);
	}

	Point(void);

};

void neighLink_xyzclcpcs(std::vector<Point> &allPoints);

void printAllNeighs(const std::vector<Point> &allPoints);

#endif // POINT_H
