#include <functional>
#include <math.h>

#include "optimise.h"

namespace thames::util::optimise{

    template<class real>
    real golden_section_search(std::function<real (real)> func, real a, real b, real tol){
        // Set parameters
        const real gr = 0.5*(sqrt(5.0) + 1);

        // Calculate test points
        real c = b - (b - a)/gr;
        real d = a + (b - a)/gr;

        // Iterate until convergence
        while (abs(b - a) > tol) {
            // Choose test points based on function values
            if (func(c) < func(d)) {
                b = d;
            } else {
                a = c;
            }

            // Recompute test points
            c = b - (b - a)/gr;
            d = a + (b - a)/gr;           
        }

        // Return converged value
        return (0.5*(a + b));
    }
    template double golden_section_search<double>(std::function<double (double)>, double, double, double);

}