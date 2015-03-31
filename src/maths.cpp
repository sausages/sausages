#include "maths.h"
#include <iostream>

std::ostream &operator<<( std::ostream &output, const vector3d &v ) { 
//	output << '(' << v.x << ',' << v.y << ',' << v.z << ')';
	output << v.x << ' ' << v.y << ' ' << v.z ;
	return output;            
}

// vector+vector
vector3d norm(const vector3d& u){
	vector3d result;
    double length;
    length = sqrt(u.x*u.x+u.y*u.y+u.z*u.z);
	result.x = u.x / length;
	result.y = u.y / length;
	result.z = u.z / length;
	return result;
}

// vector+vector
vector3d operator+(const vector3d& u, const vector3d& v){
	vector3d result;
	result.x = u.x + v.x;
	result.y = u.y + v.y;
	result.z = u.z + v.z;
	return result;
}

// vector-vector
vector3d operator-(const vector3d& u, const vector3d& v){
	vector3d result;
	result.x = u.x - v.x;
	result.y = u.y - v.y;
	result.z = u.z - v.z;
	return result;
}

// scalar*vector
vector3d operator*(const double& a, const vector3d& v){
	vector3d result;
	result.x = a * v.x;
	result.y = a * v.y;
	result.z = a * v.z;
	return result;
}

// vector/scalar
vector3d operator/(const vector3d& v, const double& a){
	vector3d result;
	result.x = v.x / a;
	result.y = v.y / a;
	result.z = v.z / a;
	return result;
}

double dot(const vector3d& u, const vector3d& v) {
	return u.x*v.x + u.y*v.y + u.z*v.z;
}

double dist(const vector3d& u, const vector3d& v) {
	double dist;
    dist = sqrt(pow(u.x-v.x,2) + pow(u.y-v.y,2) + pow(u.z-v.z,2));
    return dist;
}

double magnitude(const vector3d& u) {
	double length;
    length = sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
    return length;
}

