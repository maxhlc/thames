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

    double newton_raphson(const std::function<double (double)> &func, const std::function<double (double)> &dfunc, double xn){
        // Declare approximation variable
        double xn1;

        // Set converged flag to false
        bool converged = false;

        // Set tolerance
        double tol = 1e-10;

        // Iterate until converged
        while(!converged){
            // Update approximation
            xn1 = xn - func(xn)/dfunc(xn);

            // Converged if update is smaller than tolerance
            if(abs(xn1 - xn) < tol)
                converged = true;

            // Update previous approximation
            xn = xn1;
        }

        // Return root
        return xn1;
    }

}