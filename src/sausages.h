#ifndef FIND_SAUSAGES_H
#define FIND_SAUSAGES_H

#include <vector>
#include "point.h"

class Sausage {
	double centre_of_mass[3]; ///< Mean x/y/z of all points in the sausage
	double plane_of_best_fit[3]; ///< Unit-vector normal to plane-of-least-squares
	std::vector<double> slice_x; ///< x coord of centre of masses for slices along the sausage
	std::vector<double> slice_y; ///< y coord of centre of masses for slices along the sausage
	std::vector<double> slice_z; ///< z coord of centre of masses for slices along the sausage
	std::vector<int>slice_counter; ///< Holds # of points in slice

	public:

	Sausage(int sausageID);
	int sausageID; ///< can only be >2, doesn't make sense to have non-sausage/unsorted
	std::vector<Point> points; ///< pointers to points inside the sausage
	bool is_significant; ///< Whether the sausage is larger than some minimum size
	void find_com(); ///< Find and set centre_of_mass
	void find_pobf(); ///< Find and set plane_of_best_fit
	void estimate_sausage_length(); ///<Approximate length of sausage
	void shift_com_to_origin(std::vector<double>& xx,std::vector<double>& yy,std::vector<double>& zz); ///<Shift coordinates so that COM goes to origin
	void rotate_coord(std::vector<double>& xx,std::vector<double>& yy,std::vector<double>& zz); ///<Rotate coords to align z along direction of pobf
};

int threshold(std::vector<Point> &allPoints);
void flood_fill_separate(std::vector<Point> &allPoints, std::vector<Sausage> &allSausages);

#endif // FIND_SAUSAGES_H
