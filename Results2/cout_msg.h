#ifndef COUT_MSG_h
#define COUT_MSG_h

#include "globals.h"

// Define global parameters for messages
#define ERR  0
#define WARN  1
#define INFO  2
#define CNTR  3
#define DEB  4


#define MSG(verbosity, str)	\
	if(verbosity==ERR && verbosity<=VERBOSITY_LEVEL)   { std::cout << "\033[1;31m" << " [ERR]  :: "  << str <<  "\033[0m" << std::endl;	}	\
	if(verbosity==WARN && verbosity<=VERBOSITY_LEVEL)  { std::cout << "\033[33m"   << " [WARN] :: "  << str <<  "\033[0m" << std::endl;	}	\
	if(verbosity==INFO && verbosity<=VERBOSITY_LEVEL)  { std::cout << "\033[32m"   << " [info] :: "  << str <<  "\033[0m" << std::endl;	}	\
	if(verbosity==CNTR && verbosity<=VERBOSITY_LEVEL)  { std::cout << "\033[36m"   << " [cntr] :: "  << str <<  "\033[0m" << std::endl;	}	\
	if(verbosity==DEB && verbosity<=VERBOSITY_LEVEL)   { std::cout << "\033[94m"   << " [deb]  :: " << str <<  "\033[0m" << std::endl;	}	

#endif
