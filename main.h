#ifndef MAIN_H
#define MAIN_H

struct Point
{
	double x;  ///< cartesian x-coord of point
	double y;  ///< cartesian y-coord of point
	double z;  ///< cartesian z-coord of point
	double cl; ///< raw cl value of point
	double cp; ///< raw cp value of point
	double cs; ///< raw cl value of point
	int sausageID;
	/**< 0 for cl<threshold, 1 for unsorted.
	 * When sorted, each distinct sausage has a different number >=2.
	 * Each point is assigned to one of these sausages.
	 */
};

#endif // MAIN_H
