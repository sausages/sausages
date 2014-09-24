#include <iostream>
#include <vector>
#include <numeric>
#include <stdexcept>
#include "Eigen/Dense"
#include "io.h"
#include "point.h"
#include "sausages.h"
#include "params.h"

using namespace std;

Sausage::Sausage(int ID){
	sausageID=ID;
}

/** Sets the sausageID of each point, returning how many were below threshold
 * Set to 1 if the cl value is below the threshold (in a sausage),
 * 0 otherwise (not in a sausage)
 */
int threshold(vector<Point> &allPoints){
	int num_below_threshold=0;
	vector<Point>::iterator it;
	for (it=allPoints.begin(); it!=allPoints.end(); ++it){
		if (it->cl < params::threshold_high){
			it->sausageID=1;
			num_below_threshold++;
		} else {
			it->sausageID=0;
		}
	}
	return num_below_threshold;
}

/** Ordinary Flood-fill algorithm
 * 1) Find the first point in an unnumbered sausage
 * 2) Add it to a numbered sausage, and add all its unsorted neighbours to a stack
 * 3) While there are points in the stack, take the last one, add it to
 *     the sausage and its unsorted neighbours to the stack
 * 4) Eventually you'll have a continuous sausage, start with the next one.
 */
void flood_fill_separate(vector<Point> &allPoints, vector<Sausage> &allSausages){
	int newSausageID=2; // Start at 2 (0 is no-sausage and 1 is unsorted)
	vector<Point> stack; // Empty FILO stack, filled with points to be coloured

	for (size_t firstPoint=0; firstPoint<allPoints.size(); firstPoint++){
		// Pick the first sausage not yet coloured
		if (allPoints[firstPoint].sausageID==1){
			verbose() << "Starting sort of sausage #" << newSausageID << endl;
			allSausages.push_back(Sausage(newSausageID));
			stack.push_back(allPoints[firstPoint]);
			while (stack.size()>0){
				Point curr = stack.back();
				stack.pop_back();
				// There's got to be a nice way of doing this
				curr.self->sausageID=newSausageID;
				allSausages.back().points.push_back(curr);
				if (curr.left    && curr.left->sausageID==1)    stack.push_back(*(curr.left)) ;
				if (curr.right   && curr.right->sausageID==1)   stack.push_back(*(curr.right)) ;
				if (curr.up      && curr.up->sausageID==1)      stack.push_back(*(curr.up)) ;
				if (curr.down    && curr.down->sausageID==1)    stack.push_back(*(curr.down)) ;
				if (curr.forward && curr.forward->sausageID==1) stack.push_back(*(curr.forward)) ;
				if (curr.back    && curr.back->sausageID==1)    stack.push_back(*(curr.back)) ;
			}
			newSausageID++;
		}
	}
}

/** Classifies a given sausage
 * Four points on the sausage are found, two on either side of each colloid.
 * From each point we flood-fill towards the centre, and see which of the other points we reach
 *
 * Class 1: If we reach only the other colloid's point of the same side (i.e. above-1 -> above-2) then we
 *  have a simple ring with no twists, such as a theta configuration.
 *  Also classified as class-1 will be 2nd-loop systems where the 2nd loop is touching the 1st on only one side.
 *
 * Class 2: If we end up only on the other side of the opposite colloid (above-1 -> below-2)
 *  then we have a figure-of-eight system.
 *
 * Class 3: If we reach more than one other point we have junctions on both sides of the sausage, and are
 *  insufficiently resolved.
 *
 */
void Sausage::flood_fill_classify(void){
	/* We need to find the vector parallel to the sausage's PoBF which is perpendicular to the line joining the two colloids.
	 * Then we follow this vector out from a colloid to find suitable points on the sausage.
	 * If we arbitrarity set the vector to be length 1 we can determine it from the colloid positions and the PoBF.
	 */
	double v[3];
	double AB[3];
	double alpha,beta;

	// AB is vector joining colloids
	AB[0]=params::colloids[0][0]-params::colloids[1][0];
	AB[1]=params::colloids[0][1]-params::colloids[1][1];
	AB[2]=params::colloids[0][2]-params::colloids[1][2];
	double norm_AB=sqrt(inner_product(AB,AB+3,AB,0.0));
	debug()<<"Vector AB is {"<<AB[0]<<","<<AB[1]<<","<<AB[2]<<"} of length "<<norm_AB<<endl;

	// Plane of best fit is of form alpha*x+beta*y=z
	alpha=-plane_of_best_fit[0];
	beta=-plane_of_best_fit[1];

	// From being parallel to PoBF, perpendicular to AB and unit-vector, we can determine:
	/*
	double tmp;
	tmp  = pow(AB[1]+beta*AB[2],2) / pow(AB[0]+alpha*AB[2],2); // NB possible div0
	tmp *= (1+alpha)/(1+beta);
	v[0]=tmp/(1+tmp);
	v[1]=sqrt( (1-v[0]*v[0]*(1+alpha)) / (1+beta) );
	v[2]=alpha*v[0] + beta*v[1];
	tmp=(1+beta*beta)*( (AB[0]+alpha*AB[1])/(AB[1]+beta*AB[2]) )-alpha*beta;
	v[0]=sqrt( ( tmp*(1-beta*beta) ) / (tmp*tmp - (1+alpha*alpha+beta*beta)) );
	tmp=sqrt(v[0]*v[0]*(1+alpha*alpha+beta*beta) - beta*beta - 1 );
	v[1]= (tmp - alpha*beta*v[0]) / (1+beta*beta);
	v[2]=alpha*v[0] + beta*v[1];
	*/
	// But it's much easier to not bother normalising...
	if (abs(AB[1]+beta*AB[2]) > params::epsilon){
		v[0]=1;
		v[1]=(-alpha-AB[0])/(AB[1]+beta*AB[2]);
	}else{
		v[0]=(-AB[1]-beta*AB[2])/(AB[0]+alpha); // = 0, more or less
		v[1]=1.0/sqrt(1+beta*beta); // This pretty much normalises for free
	}
	v[2]=alpha*v[0]+beta*v[1];

	double norm_v=sqrt(inner_product(v,v+3,v,0.0));
	debug()<<"Vector v (perp. AB & || PoBF & unit) is {"<<v[0]<<","<<v[1]<<","<<v[2]<<"} of length "<<norm_v<<endl;

	// For each point in the sausage, is it between two parallel planes extended from beside a colloid, perpendicular to the line between the colloids?
	// If so, add it to a region dependant on whether it is v-wards or anti-v-wards
	vector<Point> above[2], below[2];
	for (vector<Point>::iterator it=points.begin(); it!=points.end(); it++){
		for (int iColl=0; iColl<2; iColl++){
			double p[3]; // xyz of this point
			p[0]=it->x; p[1]=it->y; p[2]=it->z;

			double u[3]; // Vector from point to colloid
			u[0]=it->x - params::colloids[iColl][0];
			u[1]=it->y - params::colloids[iColl][1];
			u[2]=it->z - params::colloids[iColl][2];

			//debug()<<"Vector u is {"<<u[0]<<","<<u[1]<<","<<u[2]<<"} of length "<<norm_u<<endl;
			//debug() << "projection: "<<projection << endl;

			/* ***********************
			 * This sweeps out a cylinder surrounding the line perpendicular to AB in the PoBF, and adds points in the cylinder to the regions.
			 * However, this may miss parts (or all) of the sausage if there are significant kinks. Better to select all points between two planes

			double distance = sqrt(norm_u*norm_u - projection*projection);
			if (distance<params::flood_fill_classify_slice_size/2){
				if (projection>0){
					// add to P1
				} else {
					// add to P2
				}
			}
			************************/

			// Eq. of plane perpendicular to AB going through point p is (x,y,z).AB - p.AB
			// Distance from point q to this plane P is |Px*qx + Py*qy + Pz*qz + P0|/norm(Px,Py,Pz)
			double first_plane[3], second_plane[3];
			for (int i=0; i<3; i++){
				first_plane[i]  = params::colloids[iColl][i] + 0.5*params::flood_fill_classify_slice_size*params::pixel_size*(AB[i]/norm_AB);
				second_plane[i] = params::colloids[iColl][i] - 0.5*params::flood_fill_classify_slice_size*params::pixel_size*(AB[i]/norm_AB);
			}
			double first_distance  = abs( inner_product(AB, AB+3, p, 0.0) - inner_product(first_plane,first_plane+3,AB,0.0) ) / norm_AB;
			double second_distance = abs( inner_product(AB, AB+3, p, 0.0) - inner_product(second_plane,second_plane+3,AB,0.0) ) / norm_AB;
			//debug()<<"distances: "<<first_distance<<","<<second_distance<<endl;

			// Is it between the two planes? If so, it is <= flood_fill_classify_slice_size from each plane
			if (first_distance < params::flood_fill_classify_slice_size*params::pixel_size &&
			   second_distance < params::flood_fill_classify_slice_size*params::pixel_size ){
				double projection = inner_product(u,u+3,v,0.0);
				if (projection>0){
					above[iColl].push_back(*(it->self));
					debug()<<it->x<<","<<it->y<<","<<it->z<<" added to region above colloid "<<iColl<<endl;
				} else {
					below[iColl].push_back(*(it->self));
					debug()<<it->x<<","<<it->y<<","<<it->z<<" added to region below colloid "<<iColl<<endl;
				}
			}
		}
	}

	verbose()<<above[0].size()<<" pixels in region above first colloid"<<endl;
	verbose()<<above[1].size()<<" pixels in region above second colloid"<<endl;
	verbose()<<below[0].size()<<" pixels in region below first colloid"<<endl;
	verbose()<<below[1].size()<<" pixels in region below second colloid"<<endl;

	/* Debug: print out all points in sausage, and points in regions */
	/*
	cout << "XXXX" << endl;
	for (vector<Point>::iterator it=points.begin(); it!=points.end(); it++){ cout << 0 <<","<<it->x<<","<<it->y<<","<<it->z<<endl; }
	for (vector<Point>::iterator it=above[0].begin(); it!=above[0].end(); it++){ cout << 1 <<","<<it->x<<","<<it->y<<","<<it->z<<endl; }
	for (vector<Point>::iterator it=below[0].begin(); it!=below[0].end(); it++){ cout << 2 <<","<<it->x<<","<<it->y<<","<<it->z<<endl; }
	for (vector<Point>::iterator it=above[1].begin(); it!=above[1].end(); it++){ cout << 3 <<","<<it->x<<","<<it->y<<","<<it->z<<endl; }
	for (vector<Point>::iterator it=below[1].begin(); it!=below[1].end(); it++){ cout << 4 <<","<<it->x<<","<<it->y<<","<<it->z<<endl; }
	cout << "XXXX" << endl;
	*/


}

void Sausage::find_com(){
	double x=0,y=0,z=0;
	vector<Point>::iterator it;
	for (it=points.begin(); it!=points.end(); ++it){
		x+=it->x;
		y+=it->y;
		z+=it->z;
	}
	centre_of_mass[0] = x/points.size();
	centre_of_mass[1] = y/points.size();
	centre_of_mass[2] = z/points.size();
	info() <<  "Centre of mass " << centre_of_mass[0] << " " << centre_of_mass[1] << " " << centre_of_mass[2] << endl;
}



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
	MatrixXd X(points.size(), 2);  ///< {x,y}-values of sausage's points, relative to sausage's CoM
	VectorXd z(points.size());      ///< z-values of sausage's points, relative to sausage's CoM
	Vector2d beta;                 ///< plane parameters

	// Centre the points around (0,0,0) (the vector normal to the plane is unaffected)
	for (size_t iPoint=0; iPoint<points.size(); iPoint++){
		X(iPoint,0) = points[iPoint].x - centre_of_mass[0];
		X(iPoint,1) = points[iPoint].y - centre_of_mass[1];
		z(iPoint) = points[iPoint].z - centre_of_mass[2];
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

void Sausage::rotate_to_xy_plane(double** pointsArray){

	debug() << "Rotating to xy-plane" << endl;

	// loop over all coordinate points and rotate
	for(size_t i = 0; i != points.size(); i++) {
		pointsArray[i][0] = rotation_matrix(0)*pointsArray[i][0] + rotation_matrix(1)*pointsArray[i][1] + rotation_matrix(2)*pointsArray[i][2];
		pointsArray[i][1] = rotation_matrix(3)*pointsArray[i][0] + rotation_matrix(4)*pointsArray[i][1] + rotation_matrix(5)*pointsArray[i][2];
		pointsArray[i][2] = rotation_matrix(6)*pointsArray[i][0] + rotation_matrix(7)*pointsArray[i][1] + rotation_matrix(8)*pointsArray[i][2];
	}
	return;
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

void Sausage::calculate_sausage_length(double **slice_positions){

	// loop over all COM of slices and calculate length
	length = 0;
	for(int i = 0; i != nSlices-1; i++) {
		length += sqrt( pow(slice_positions[i][0]-slice_positions[i+1][0],2)+
				pow(slice_positions[i][1]-slice_positions[i+1][1],2)+
				pow(slice_positions[i][2]-slice_positions[i+1][2],2));
	}
	length += sqrt( pow(slice_positions[nSlices-1][0]-slice_positions[0][0],2)+
			pow(slice_positions[nSlices-1][1]-slice_positions[0][1],2)+
			pow(slice_positions[nSlices-1][2]-slice_positions[0][2],2));

	info() << "The estimated length of the sausage is " << length << endl;
	return;
}

void Sausage::estimate_sausage_length(){

	info() << "--------> Estimating sausage length... " << endl;

	// Make array of points, shifted so that COM is in origin
	double **rotatedPoints = new double*[points.size()];
	for (size_t i=0; i<points.size(); i++){
		rotatedPoints[i] = new double[3];
		rotatedPoints[i][0]=points[i].x - centre_of_mass[0];
		rotatedPoints[i][1]=points[i].y - centre_of_mass[1];
		rotatedPoints[i][2]=points[i].z - centre_of_mass[2];
	}

	// calculate rotation matrix and its inverse, so that plane of best fit projects onto xy-plane
	calculate_rotation_matrix();
	rotate_to_xy_plane(rotatedPoints);

	// Find number of slices
	nSlices = points.size()/params::points_per_slice;
	info() << "Sausage was divided in " << nSlices << " slices." << endl;

	double **slice_positions = new double*[nSlices];
	for (int i=0; i<nSlices; i++){
		slice_positions[i] = new double[3];
	}
	double *slice_counter = new double[nSlices];


	// find what slice we are in (loop over all points)
	double theta;
	int sliceId;
	for(size_t i = 0; i != points.size(); i++) {

		theta=atan2(rotatedPoints[i][1],rotatedPoints[i][0])+M_PI;
		sliceId = (int) nSlices*theta/(2.0*M_PI);
		if (sliceId < 0 || sliceId >= nSlices ) {
			cerr << "SliceId is out of bound" << endl;
			exit(EXIT_FAILURE);
		}

		//add x,y,z coordinates to that slice
		slice_positions[sliceId][0]+=rotatedPoints[i][0];
		slice_positions[sliceId][1]+=rotatedPoints[i][1];
		slice_positions[sliceId][2]+=rotatedPoints[i][2];
		slice_counter[sliceId]++;
	}

	// work out centre of mass for each slice
	for(int k = 0; k != nSlices; k++) {
		if (slice_counter[k] == 0) {
			cerr << "Slice is empty!!" << endl;
			exit(EXIT_FAILURE);}
		else{
			slice_positions[k][0] = slice_positions[k][0]/slice_counter[k];
			slice_positions[k][1] = slice_positions[k][1]/slice_counter[k];
			slice_positions[k][2] = slice_positions[k][2]/slice_counter[k];

			debug() << slice_counter[k] << " COM " << slice_positions[k][0] << " " << slice_positions[k][1] << " " << slice_positions[k][2] << endl;}
	}

	// convert COM's of slice back to initial frame
	//rotate_from_xy_plane(slice_positions);

	// calculate length
	calculate_sausage_length(slice_positions);

	// clean up
	debug() << "deleting rotated_points" << endl;
	for (size_t i=0; i<points.size(); i++){
		delete rotatedPoints[i];
	}
	delete rotatedPoints;

	debug() << "deleting slice_positions" << endl;
	for (int i=0; i<nSlices; i++){
		delete slice_positions[i];
	}
	delete slice_positions;

	return;
}
