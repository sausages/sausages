#ifndef MATHS_H
#define MATHS_H

#include <iostream>
#include <math.h>


typedef struct vector3d_ {
	double x,y,z;
} vector3d;

// <<normalised vector
vector3d norm(const vector3d& u);

// <<vector
std::ostream &operator<<( std::ostream &output, const vector3d &v );

// vector+vector
vector3d operator+(const vector3d& u, const vector3d& v);

// vector-vector
vector3d operator-(const vector3d& u, const vector3d& v);

// scalar*vector
vector3d operator*(const double& a, const vector3d& v);

// vector/scalar
vector3d operator/(const vector3d& v, const double& a);

double dot(const vector3d& u, const vector3d& v);

double dist(const vector3d& u, const vector3d& v);

double magnitude(const vector3d& u);

#endif // MATHS_H
