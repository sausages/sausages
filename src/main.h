#ifndef MAIN_H
#define MAIN_H

#include <vector>
#include "maths.h"

class Point;
class Sausage;

class Model {
	public:
	int total_points; ///< How many points are in the input file?
	std::vector<int> threshold_points; ///< Which points are below the threshold? Stores their allPointsIndex
	std::vector<Vector3d> colloidPos; ///< Positions of the colloids
	std::vector<double> colloidRadii; ///< Radii of the colloids
	std::vector<Point*> allPoints; ///< All points in simulation
	std::vector<Sausage> allSausages; ///< All sausages in simulation
};


// There doesn't seem to be an existing neat way of pythonic "element in vec"
// So here is something that will allow "elementOf(vec,element)"
template <class T>
bool elementOf(std::vector<T> vec , T element){
	return (std::find(vec.begin(), vec.end(), element) != vec.end());
}



#endif // MAIN_H
