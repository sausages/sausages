#ifndef MAIN_H
#define MAIN_H

// enums are implicitly assigned to 0,1,2.. so e.g. DEBUG > WARNING
enum verbosityLevel {ERROR, WARNING, NORMAL, VERBOSE, DEBUG};

namespace params{
	extern const verbosityLevel verbosity;
	extern const double threshold_level;
	extern const double silent_ignore_size;
	extern const double min_sausage_size;
}


#endif // MAIN_H
