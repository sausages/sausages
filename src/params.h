#ifndef PARAMS_H
#define PARAMS_H

// enums are implicitly assigned to 0,1,2.. so e.g. DEBUG > WARNING
enum verbosityLevel {ERROR, WARNING, INFO, VERBOSE, DEBUG};

/*
namespace params{
	extern const verbosityLevel verbosity;
	extern const double threshold_high; ///< Sausages should be without gaps at this level
	extern const double threshold_low; ///< All different defect structures should be distinguishable (i.e. no ambigious blobs)
	extern const double silent_ignore_size;
	extern const double min_sausage_size;
	extern const int points_per_slice;
}
*/
namespace params{
	extern verbosityLevel verbosity;
	extern double threshold_high; ///< Sausages should be without gaps at this level
	extern double threshold_low; ///< All different defect structures should be distinguishable (i.e. no ambigious blobs)
	extern double silent_ignore_size;
	extern double min_sausage_size;
	extern int points_per_slice;
}

void set_params(char *paramFile);

#endif // PARAMS_H
