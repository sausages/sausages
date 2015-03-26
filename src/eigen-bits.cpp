#include "Eigen/Dense"
#include "io.h"
#include "sausages.h"

using namespace std;

void Sausage::find_pobf(){
	/**
	 * If we treat x,y,z as the distance of a point from the centre-of-mass of its sausage, then:
	 *
	 * \f[
	 * \vec{x_i} \cdot \vec{\beta} = z_i
	 * \bf{X} \cdot \vec{\beta} = \vec{z}
	 * \f]
	 *
	 * The optimal fit \f$\vec{\hat{\beta}}\f$ is:
	 * \f[
	 * \vec{\hat{\beta}} = (\bf{X}^T \bf{X} )^{-1} \bf{X}^T \vec{z}
	 * \f]
	 */


	using namespace Eigen;

	if (have_pobf) return;
	have_pobf = true;

	MatrixXd X(points.size(), 2);  ///< {x,y}-values of sausage's points, relative to sausage's CoM
	VectorXd z(points.size());      ///< z-values of sausage's points, relative to sausage's CoM
	Vector2d beta;                 ///< plane parameters

	// Centre the points around (0,0,0) (the vector normal to the plane is unaffected)
	for (size_t iPoint=0; iPoint<points.size(); iPoint++){
		X(iPoint,0) = points[iPoint]->x - centre_of_mass[0];
		X(iPoint,1) = points[iPoint]->y - centre_of_mass[1];
		z(iPoint) = points[iPoint]->z - centre_of_mass[2];
	}

	// Eigen can give us a Jacobi SVD, and then use that SVD to give a (guaranteed least-squares) solution of Ax=b
	beta = X.jacobiSvd(ComputeThinU | ComputeThinV).solve(z) ;
	debug() << "Least-squares plane is: " << beta(0) << "x + " << beta(1) << "y = z" << endl;

	// Points O=(0,0,0) a=(1,0,alpha) b=(0,1,beta) are in the plane
	// Therefore (1,0,alpha) x (0,1,beta) = (-alpha, -beta, 1) is prependicular to the plane
	plane_of_best_fit[0]=-beta[0];
	plane_of_best_fit[1]=-beta[1];
	plane_of_best_fit[2]=1;
	debug() << "Vector perpendicular to the plane is (" << -beta[0] << "," << -beta[1] << ",1)" << endl;

}

/**
 * We want a rotation matrix which will rotate the sausage, once translated to the origin,
 * so that its PoBF is the xy-plane.
 * To do this, find the matrix which rotates the vector normal to the plane to lie along the z-axis.
 * Using the Rodrigues' Rotation Formula, the matrix which rotates \f$\hat{a}\f$ onto \f$\hat{b}\f$ is:
 * \f[
 * R = I + s [v]_\times + (1-c) [v]_\times^2
 * \f]
 * where \f$v=\hat{a}\times \hat{b}\f$ , \f$s=sin(\theta)=|\hat{a}\times\hat{b}|\f$ and \f$c=cos{\theta}=\hat{a}\cdot\hat{b}\f$
 * and
 * \f[
 * [v]_times =
 * \begin{pmatrix}
 *  0   & -v_3 &  v_2 \\
 *  v_3 &  0   & -v_2 \\
 * -v_2 &  v_2 &  0
 *  \end{pmatrix}
 * \f]
 * is the cross-product matrix of \f$v\f$
 *
 */
void Sausage::calculate_rotation_matrix(void){

	if (have_rotation_matrix) return;
	have_rotation_matrix = true;

	find_pobf();

	Eigen::Vector3d a(plane_of_best_fit[0],plane_of_best_fit[1],plane_of_best_fit[2]);
	Eigen::Vector3d z(0,0,1);
	a.normalize();

	Eigen::Vector3d v=a.cross(z);

	double s=v.norm();
	double c=a.dot(z);

	Eigen::Matrix3d v_cross_mat;
	v_cross_mat << 0   , -v(2),  v(1),
	        v(2),  0   , -v(0),
	       -v(1),  v(0),  0   ;

	rotation_matrix = Eigen::Matrix3d::Identity() + s*v_cross_mat + (1.0-c)*v_cross_mat*v_cross_mat;

	debug() << "Rotation matrix is:" << endl << rotation_matrix << endl << endl;

	return;
}
