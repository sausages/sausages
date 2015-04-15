#ifndef PARAMS_H
#define PARAMS_H

#include <vector>
#include "main.h"
#include "maths.h"

// enums are implicitly assigned to 0,1,2.. so e.g. DEBUG > WARNING
enum verbosityLevel {ERROR, WARNING, INFO, VERBOSE, DEBUG};

namespace params{
	extern verbosityLevel verbosity;
	extern std::string brief_filename; ///< 'brief' file is for standardised output, in a different file to cout
	extern int brief_version; ///< each brief version is standardised
	extern std::string thresholded_filename; ///< 'thresholded' file is for output of points below the threshold, in a different file to cout
	extern std::string prejoin_sausage_filename; /// file for output of points in sausages before joining, one file per sausage
	extern std::string sausage_filename; /// 'sausage' file is for output of points in sausages after joining endpoints, one file per sausage
	extern double threshold; ///< c_l threshold for includion into a sausage. All different defect structures should be distinguishable (i.e. no ambigious blobs)
	extern double silent_ignore_size; ///< If a sausage is smaller than this fraction of all points below the threshold, silently ignore it
	extern double min_sausage_size; ///< A sausage is only 'significant' if it is larger than this fraction of all points below the threshold
	extern double pixel_size; ///< Distance between nearest-neighbour points
	extern int points_per_slice; ///< How many points (on average) are in each slice of the 'pearl-necklace' sausage-length measurer
	extern bool colloidsInParamFile; ///< Colloids shouldn't be both in param and input file
	extern double flood_fill_classify_slice_size; ///< How many pixels wide should the regions in flood_fill_classify be?
	extern double ratio_two_rings; ///< If ratio of two relevant sausage is below this, it's the two ring structure
	extern double ratio_2nd_loop; ///< If the ratio of two relevant sausages is above this, it's the 2nd loop structure
	extern double epsilon; ///< Some small number for comparison to zero
	extern double max_sausage_gap_size; ///< Maximum distance (in units of pixel-length) between endpoints that will be joined by join_endpoints()
	extern double min_sausage_gap_next_neigh_distance; ///< If there are multiple neighbouring endpoints within this radius (in units of pixels), we throw an error as it's too close to call between candidates
}


void set_params(char *paramFile, std::vector<Vector3d> &colloidPos);

#endif // PARAMS_H
