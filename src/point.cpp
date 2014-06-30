#include <iostream>
#include <cmath>
#include "point.h"

using namespace std;

/** Links each point in allPoints to its 6 nearest-neighbours.
 * x:+right/-left ; y:+up/-down ; z:+forward/-back
 * @warning Assumes that points are listed in a specific order:
 * @warning   z varies first (increasing), then y, then x.
 * @warning Also assumes that grid is cubic, i.e. N=len(x)=len(y)=len(z)
 */
void neighLink_xyzclcpcs(vector<Point> &allPoints){
	unsigned int N = round(pow(allPoints.size(),1.0/3.0));///< side-length of cube
	// Don't use iterators, indices are important here
	for (size_t i=0; i<allPoints.size(); i++){
		// Pointer to self, pretty sure this is irrelevant, but
		// I currently need it for the flood-fill (iterators are copies)
		allPoints[i].self = &(allPoints[i]);
		// If not right-most
		if (((i/(N*N))+1)%N != 0){
			allPoints[i].right  = &(allPoints[i+N*N]);
			allPoints[i+N*N].left = &(allPoints[i]);
		}
		// If not upper-most
		if (((i/N)+1)%N != 0){
			allPoints[i].up   = &(allPoints[i+N]);
			allPoints[i+N].down = &(allPoints[i]);
		}
		// If not forward-most
		if ((i+1)%N != 0){
			allPoints[i].forward = &(allPoints[i+1]);
			allPoints[i+1].back   = &(allPoints[i]);
		}
	}
}


void printAllNeighs(const vector<Point> &allPoints){
	Point p;
	const size_t max_PointsToPrint=100;
	size_t num_PointsToPrint=min(allPoints.size(),max_PointsToPrint);
	for (size_t i = 0; i<num_PointsToPrint ; i++){
		p=allPoints[i];
		cout << i << endl;
		cout << "I am " << &allPoints[i] << endl; // Don't print &p
		cout << "at x:" << p.x << " y:" << p.y << " z:" << p.z << endl;
		cout << "l " << p.left << endl;
		cout << "r " << p.right << endl;
		cout << "u " << p.up   << endl;
		cout << "d " << p.down << endl;
		cout << "f " << p.forward << endl;
		cout << "b " << p.back << endl;
		cout << endl;
	}
}
