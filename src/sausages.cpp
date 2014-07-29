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

	//plane_of_best_fit

	info() << "Least-squares plane is " << X.jacobiSvd(ComputeThinU | ComputeThinV).solve(z) << endl;

}

void Sausage::shift_com_to_origin(vector<double>& xx, vector<double>& yy, vector<double>& zz){

	int l=0;
	vector<Point>::iterator it;
	for (it=points.begin(); it!=points.end(); ++it){
		xx[l] = it->x - centre_of_mass[0]; 
		yy[l] = it->y - centre_of_mass[1]; 
		zz[l] = it->z - centre_of_mass[2]; 
		l++;
	}
	return;
}
void Sausage::rotate_coord(vector<double>& xx, vector<double>& yy, vector<double>& zz, double A[9]){

	// loop over all coordinate points and rotate
	for(std::vector<int>::size_type i = 0; i != xx.size(); i++) {
	   	xx[i] = A[0]*xx[i] + A[1]*yy[i] + A[2]*zz[i];
	   	yy[i] = A[3]*xx[i] + A[4]*yy[i] + A[5]*zz[i];
	   	zz[i] = A[6]*xx[i] + A[7]*yy[i] + A[8]*zz[i];
	    }
	return;
}

void Sausage::calculate_rotation_matrix(double alpha, double beta){

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

void Sausage::estimate_sausage_length(){

	info() << "--------> Estimating sausage length... " << endl;
 	int nop = points.size(); // number of points in sausage
	std::vector<double>xx(nop); // arrays for temporary coordinates
	std::vector<double>yy(nop);
	std::vector<double>zz(nop);
	// find centre of mass
	find_com();
	// shift coordinates so that COM is in origin
	shift_com_to_origin(xx,yy,zz);
	
	/*for(std::vector<int>::size_type i = 0; i != xx.size(); i++) {
		info() << xx[i] << " " << yy[i] << " " << zz[i] << endl;
	}*/
	
	// find plane of best fit
	find_pobf();
	alpha =	0.124972; 
	beta = 0.133109;
	// calculate rotation matrix and it's inverse, so that plane of best fit projects onto xy-plane
	calculate_rotation_matrix(alpha,beta);
	rotate_coord(xx,yy,zz,rotation_matrix);
	// Find number of slices
	int nos = (int)points.size()/100;
	info() << "Sausage was divided in " << nos << " slices." << endl;
	slice_counter.resize(nos);
	slice_x.resize(nos);slice_y.resize(nos);slice_z.resize(nos);
	// find what slice we are in (loop over all points)
	double theta; 
	int sliceId;
	for(std::vector<int>::size_type i = 0; i != xx.size(); i++) {
		// Is there fn for pi or does it have to be defined as a constant?
		theta=atan2(yy[i],xx[i])+3.14159265358;
		sliceId = (int) (double)nos*theta/(2.0*3.14159265358);
		if (sliceId < 0 || sliceId >= nos ) {
			cerr << "SliceId is out of bound" << endl;
			exit(EXIT_FAILURE);}
		//add x,y,z coordinates to that slice
		slice_x[sliceId]+=xx[i]; 
		slice_y[sliceId]+=yy[i];
		slice_z[sliceId]+=zz[i];
		slice_counter[sliceId]++;
	}
	// work out centre of mass for each slice	
	for(std::vector<int>::size_type k = 0; k != slice_x.size(); k++) {
		if (slice_counter[k] == 0) {
			cerr << "Slice is empty!!" << endl;
			exit(EXIT_FAILURE);}
		else{
			slice_x[k] = slice_x[k]/slice_counter[k];
			slice_y[k] = slice_y[k]/slice_counter[k];
			slice_z[k] = slice_z[k]/slice_counter[k];
			// print centre of mass for all slices
			debug() << slice_counter[k] << " COM " << slice_x[k] << " " << slice_y[k] << " " << slice_z[k] << endl;}
	}
	// convert COM's of slice back to initial frame
	rotate_coord(xx,yy,zz,inv_rotation_matrix);
	return;
}

