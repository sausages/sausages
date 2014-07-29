#include <iostream>
#include <vector>
#include <stdexcept>
#include "Eigen/Dense"
#include "io.h"
#include "point.h"
#include "sausages.h"
#include "main.h"

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
		if (it->cl < params::threshold_level){
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
	 * \vec{x_i} \dot \vec{\beta} = z_i
	 * \bf{X} \dot \vec{\beta} = \vec{z}
	 *
	 * The optimal fit \vec{\hat{\beta}} is:
	 * \vec{\hat{\beta}} = (\bf{X}^T \bf{X} )^{-1} \bf{X}^T \vec{z}
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
	debug() << "Least-squares plane is " << endl << beta << endl;

	// Points O=(0,0,0) a=(1,0,alpha) b=(0,1,beta) are in the plane
	// Therefore (1,0,alpha) x (0,1,beta) = (-alpha, -beta, 1) is prependicular to the plane
	plane_of_best_fit[0]=-beta[0];
	plane_of_best_fit[1]=-beta[1];
	plane_of_best_fit[2]=1;
	debug() << "Vector perpendicular to the plane is (" << -beta[0] << "," << -beta[1] << ",1)" << endl;

}

void Sausage::rotate_to_xy_plane(double** pointsArray){

	// loop over all coordinate points and rotate
	for(size_t i = 0; i != points.size(); i++) {
		pointsArray[i][0] = rotation_matrix[0]*pointsArray[i][0] + rotation_matrix[1]*pointsArray[i][1] + rotation_matrix[2]*pointsArray[i][2];
		pointsArray[i][1] = rotation_matrix[3]*pointsArray[i][0] + rotation_matrix[4]*pointsArray[i][1] + rotation_matrix[5]*pointsArray[i][2];
		pointsArray[i][1] = rotation_matrix[6]*pointsArray[i][0] + rotation_matrix[7]*pointsArray[i][1] + rotation_matrix[8]*pointsArray[i][2];
	}
	return;
}

void Sausage::calculate_rotation_matrix(void){

	// loop over all coordinate points and rotate
/*	for(std::vector<int>::size_type i = 0; i != xx.size(); i++) {
		xx[i] = plane_of_best_fit[0]*zz[i];
		yy[i] = plane_of_best_fit[1]*zz[i];
		zz[i] = plane_of_best_fit[2]*zz[i];
	}*/
	//dummy matrix
	rotation_matrix[0] = 1.0; rotation_matrix[1] = 0.0; rotation_matrix[2] = 0.0;
	rotation_matrix[3] = 0.0; rotation_matrix[4] = 1.0; rotation_matrix[5] = 0.0;
	rotation_matrix[6] = 0.0; rotation_matrix[7] = 0.0; rotation_matrix[8] = 1.0;

	inv_rotation_matrix[0] = 1.0; inv_rotation_matrix[1] = 0.0; inv_rotation_matrix[2] = 0.0;
	inv_rotation_matrix[3] = 0.0; inv_rotation_matrix[4] = 1.0; inv_rotation_matrix[5] = 0.0;
	inv_rotation_matrix[6] = 0.0; inv_rotation_matrix[7] = 0.0; inv_rotation_matrix[8] = 1.0;

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
	nSlices = points.size()/params::pointsPerSlice;
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

	return;
}
