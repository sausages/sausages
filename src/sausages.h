#ifndef FIND_SAUSAGES_H
#define FIND_SAUSAGES_H

#include <vector>
#include "Eigen/Dense"
#include "point.h"

// All Eigen-related methods are defined in eigen-bits.cpp

class Sausage {
	double centre_of_mass[3]; ///< Mean x/y/z of all points in the sausage

	double **slice_positions; ///< xyz positions of COM of points within a slice
	int *slice_counter; ///< Holds # of points in slice
	int nSlices; ///< # slices
	double length; ///< Estimated length of sausage

	Eigen::Matrix3d rotation_matrix; //rotation matrix from initial to new frame

	void rotate_to_xy_plane(double** points); ///<Rotate coords (rotation matrix*coord)
	void rotate_from_xy_plane(double** points); ///<Rotate coords (rotation matrix*coord)
	void calculate_rotation_matrix(void); ///<Calculates rotation matrix and its inverse
	void calculate_sausage_length(double** slice_positions); ///<Calculate the length using COMs of slices

	public:

	Sausage(int sausageID); ///< Minimal constructor

	int sausageID; ///< can only be >2, doesn't make sense to have non-sausage/unsorted
	std::vector<Point> points; ///< pointers to points inside the sausage
	bool is_significant; ///< Whether the sausage is larger than some minimum size

	void find_com(void); ///< Find and set centre_of_mass
	void find_pobf(void); ///< Find and set plane_of_best_fit
	void estimate_sausage_length(void); ///<Approximate length of sausage
	void flood_fill_classify(void);
	void find_endpoints(void);
	double plane_of_best_fit[3]; ///< Unit-vector normal to plane-of-least-squares
};


/* Functions below are not particular to a specific sausage */
int threshold(std::vector<Point> &allPoints); ///< Pick out all sausages with cl>threshold
void flood_fill_separate(std::vector<Point> &allPoints, std::vector<Sausage> &allSausages); ///< use flood-fill to distinguish between non-touching sausages

#endif // FIND_SAUSAGES_H
