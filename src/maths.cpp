#include "maths.h"
#include <cmath>
#include <iostream>

std::ostream &operator<<( std::ostream &output, const Vector3d &v ) {
	output << '(' << v.x << ',' << v.y << ',' << v.z << ')';
	return output;            
}

// vector+vector
Vector3d operator+(const Vector3d& u, const Vector3d& v){
	return Vector3d (u.x+v.x, u.y+v.y, u.z+v.z);
}

// vector-vector
Vector3d operator-(const Vector3d& u, const Vector3d& v){
	return Vector3d (u.x-v.x, u.y-v.y, u.z-v.z);
}

// scalar*vector
Vector3d operator*(const double& a, const Vector3d& v){
	return Vector3d (a*v.x, a*v.y, a*v.z);
}

// vector/scalar
Vector3d operator/(const Vector3d& v, const double& a){
	return Vector3d (v.x/a, v.y/a, v.z/a);
}

double dot(const Vector3d& u, const Vector3d& v) {
	return u.x*v.x + u.y*v.y + u.z*v.z;
}

double mag(const Vector3d& u){
	return sqrt(dot(u,u));
}
