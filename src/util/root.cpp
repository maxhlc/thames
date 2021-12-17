#include <functional>
#include <math.h>

#include "root.h"
#include "optimise.h"

namespace thames::util::root{

    double golden_section_search(std::function<double (double)> func, double a, double b){
        // Declare function for minimisation (root at minimum of absolute of the function)
        std::function<double (double)> f = [func](double x) {return abs(func(x));};

        // Minimise function using the golden section search
        double x0 = thames::util::optimise::golden_section_search(f, a, b);

        // Return root
        return x0;
    }

}