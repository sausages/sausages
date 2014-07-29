#ifndef PARAMS_H
#define PARAMS_H

// enums are implicitly assigned to 0,1,2.. so e.g. DEBUG > WARNING
enum verbosityLevel {ERROR, WARNING, INFO, VERBOSE, DEBUG};

namespace params{
	extern const verbosityLevel verbosity;
	extern const double threshold_level;
	extern const double silent_ignore_size;
	extern const double min_sausage_size;
	extern const int pointsPerSlice;
}

void set_params(char *paramFile);

#endif // PARAMS_H
