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
	 * In practice we should use this, not explicit left/right etc
	 * The exception is when reading input/setting up
	 * Advantage of array is possibility of looping over it
	 */
	Point* neighbours[6];

	// Pointer to self, pretty sure this is irrelevant, but...
	// I currently need it for the flood-fill (iterators are copies) and == operator
	Point* self;

	// Assume if our self-pointer is the same then we are the same
	// Currently needed for the find algorithm in flood-fill
	bool operator == (const Point &Ref) const {
		return(this->self == Ref.self);
	}
	// Sort points (for uniqueness) by value of pointer
	bool operator < (const Point &Ref) const {
		return(this->self < Ref.self);
	}

};

void neighLink_xyzclcpcs(std::vector<Point> &allPoints);

void printAllNeighs(const std::vector<Point> &allPoints);

#endif // POINT_H
