#ifndef MAIN_H
#define MAIN_H

#include <vector>
#include "maths.h"
#include "point.h"
#include "sausages.h"

class Model {
	public:
	int total_points; ///< How many points are in the input file?
	int threshold_points; ///< How many points are below the threshold?
	std::vector<vector3d> colloidPos; ///< Positions of the colloids
	std::vector<double> colloidRadii; ///< Radii of the colloids
	std::vector<Point> allPoints; ///< All points in simulation
	std::vector<Sausage> allSausages; ///< All sausages in simulation
};

/*
namespace model{
	extern std::vector<Point> allPoints; ///< All points in simulation
	extern std::vector<Sausage> allSausages; ///< All sausages in simulation

	extern std::vector<vector3d> colloidPos; ///< Positions of the colloids
	extern std::vector<double> colloidRadii; ///< Radii of the colloids

}
*/

// There doesn't seem to be an existing neat way of pythonic "element in vec"
// So here is something that will allow "elementOf(vec,element)"
template <class T>
bool elementOf(std::vector<T> vec , T element){
	return (std::find(vec.begin(), vec.end(), element) != vec.end());
}



#endif // MAIN_H
