#ifndef MATHS_H
#define MATHS_H

#include <iostream>

class Vector3d {
	public:
	double x,y,z;
	// constructor
	Vector3d(double x, double y, double z): x(x),y(y),z(z) {};
};


// <<vector
std::ostream &operator<<( std::ostream &output, const Vector3d &v );

// vector+vector
Vector3d operator+(const Vector3d& u, const Vector3d& v);

// vector-vector
Vector3d operator-(const Vector3d& u, const Vector3d& v);

// scalar*vector
Vector3d operator*(const double& a, const Vector3d& v);

// vector/scalar
Vector3d operator/(const Vector3d& v, const double& a);

double dot(const Vector3d& u, const Vector3d& v);


#endif // MATHS_H
