#ifndef UTIL_OPTIMISE
#define UTIL_OPTIMISE

#include <functional>
#include <math.h>

namespace thames::util::optimise{

    double golden_section_search(std::function<double (double)> func, double a, double b);

}

#endif