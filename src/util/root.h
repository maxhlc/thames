#ifndef UTIL_ROOT
#define UTIL_ROOT

#include <functional>
#include <math.h>

namespace thames::util::root{

    double golden_section_search(std::function<double (double)> func, double a, double b);

}

#endif