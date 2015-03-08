#ifndef MAIN_H
#define MAIN_H

#include <vector>
#include "maths.h"
#include "point.h"
#include "sausages.h"

namespace model{
	extern std::vector<Point> allPoints; ///< All points in simulation
	extern std::vector<Sausage> allSausages; ///< All sausages in simulation

	extern std::vector<vector3d> colloidPos; ///< Positions of the colloids
	extern std::vector<double> colloidRadii; ///< Radii of the colloids

}



#endif // MAIN_H
