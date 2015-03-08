#include "maths.h"
#include <iostream>

std::ostream &operator<<( std::ostream &output, const vector3d &v ) { 
	output << '(' << v.x << ',' << v.y << ',' << v.z << ')';
	return output;            
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


