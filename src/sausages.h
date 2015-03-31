#ifndef FIND_SAUSAGES_H
#define FIND_SAUSAGES_H

#include <vector>
#include "Eigen/Dense"
#include "maths.h"
#include "point.h"

// All Eigen-related methods are defined in eigen-bits.cpp

class Sausage {
	double centre_of_mass[3]; ///< Mean x/y/z of all points in the sausage
	double length; ///< Estimated length of sausage

	bool have_pobf; ///< have we calculated the plane of best fit yet?
	bool have_rotation_matrix; ///< have we calculated the rotation matrix yet?
	Eigen::Matrix3d rotation_matrix; ///< rotation matrix from initial to new frame

	void rotate_to_xy_plane(double** points); ///<Rotate coords (rotation matrix*coord)
	void rotate_from_xy_plane(double** points); ///<Rotate coords (rotation matrix*coord)
	void calculate_rotation_matrix(void); ///<Calculates rotation matrix and its inverse
    void calculate_com_sphere(vector3d center, double radius, vector3d &com); ///< Calculate COM for 10 points nearest to a point com_x,com_y,com_z, used for halfsphere_tracking
//    void calculate_sausage_length_halfsphere(void); ///< Calculate length by connecting all COMs of halfsphere tracking by straight lines

	public:

	Sausage(int sausageID); ///< Minimal constructor

	int sausageID; ///< Which sausage this is
	std::vector<Point*> points; ///< pointers to points inside the sausage
	size_t endpoints[2]; ///< Indices of the points in self.points which correspond to the endpoints of the sausage
	bool is_significant; ///< Whether the sausage is larger than some minimum size
	std::vector<vector3d> halfsphere_COMs; ///< COMs produced by halfsphere tracking

	void find_com(void); ///< Find and set centre_of_mass
	void find_pobf(void); ///< Find and set plane_of_best_fit
    void halfsphere_tracking(); ///< Tracks sausage and saves all coms 
	void flood_fill_classify(const std::vector<vector3d> colloidPos);
	void find_endpoints(void);
	double plane_of_best_fit[3]; ///< Unit-vector normal to plane-of-least-squares
};


/* Functions below are not particular to a specific sausage */
int threshold(std::vector<Point> &allPoints); ///< Pick out all sausages with cl>threshold
void flood_fill_separate(std::vector<Point> &allPoints, std::vector<Sausage> &allSausages); ///< use flood-fill to distinguish between non-touching sausages
void join_endpoints(std::vector<Sausage> &allSausages, std::vector<int> &relevantSausages); ///< Artificially join gaps between endpoints less than params::max_endpoint_separation

#endif // FIND_SAUSAGES_H
