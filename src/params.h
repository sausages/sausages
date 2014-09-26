#ifndef PARAMS_H
#define PARAMS_H

// enums are implicitly assigned to 0,1,2.. so e.g. DEBUG > WARNING
enum verbosityLevel {ERROR, WARNING, INFO, VERBOSE, DEBUG};

namespace params{
	extern verbosityLevel verbosity;
	extern double threshold; ///< Sausages should be without gaps at this level
	extern double silent_ignore_size; ///< If a sausage is smaller than this fraction of all points below the threshold, silently ignore it
	extern double min_sausage_size; ///< A sausage is only 'significant' if it is larger than this fraction of all points below the threshold
	extern double pixel_size; ///< Distance between nearest-neighbour points
	extern int points_per_slice; ///< How many points (on average) are in each slice of the 'pearl-necklace' sausage-length measurer
	extern double colloids[2][3]; ///< xyz positions of the two colloids
	extern double flood_fill_classify_slice_size; ///< How many pixels wide should the regions in flood_fill_classify be?
	extern double ratio_two_rings; ///< If ratio of two relevant sausage is below this, it's the two ring structure
	extern double ratio_2nd_loop; ///< If the ratio of two relevant sausages is above this, it's the 2nd loop structure
	extern double epsilon; ///< Some small number for comparison to zero
}


void set_params(char *paramFile);

#endif // PARAMS_H
