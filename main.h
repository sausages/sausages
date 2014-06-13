#ifndef MAIN_H
#define MAIN_H

// enums are implicitly assigned to 0,1,2.. so e.g. DEBUG > WARNING
enum verbosityLevel {ERROR, WARNING, NORMAL, VERBOSE, DEBUG};

namespace params{
	extern verbosityLevel verbosity;
}


#endif // MAIN_H
