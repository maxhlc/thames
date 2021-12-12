#include <functional>
#include <math.h>

#include "optimise.h"

namespace thames::util::optimise{

    double golden_section_search(std::function<double (double)> func, double a, double b){
        // TODO: documentation

        // Set parameters
        const double gr = 0.5*(sqrt(5.0) + 1);
        const double tol = 1e-14;

        // Calculate test points
        double c = b - (b - a)/gr;
        double d = a + (b - a)/gr;

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

}