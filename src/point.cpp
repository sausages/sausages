#include <iostream>
#include "point.h"
#include "params.h"
#include "io.h"

using namespace std;

Point::Point(void){
	left=right=up=down=forward=back=NULL;
	sausageID=-1; // unsorted
	isInASausage=false;
}

Point::operator Vector3d(){
	Vector3d v (x,y,z);
	return v;
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
