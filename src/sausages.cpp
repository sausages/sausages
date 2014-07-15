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

void Sausage::rotate_coord(vector<double>& xx, vector<double>& yy, vector<double>& zz){

	// loop over all coordinate points and rotate
	for(std::vector<int>::size_type i = 0; i != xx.size(); i++) {
	   // 	xx[i] = ...; plane_of_best_fit[3]
	    }
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
	// find plane of best fit
	find_pobf();
	// rotate coordinate system so that plane of best fit projects onto xy-plane
	rotate_coord(xx,yy,zz);
	// Find number of slices
	int nos = (int)points.size()/100;
	info() << "Sausage was divided in " << nos << " slices." << endl;
	// find which slice we are in
	// work out centre of mass for each slice
	return;
}


