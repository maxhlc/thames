#include <functional>
#include <math.h>

#include "optimise.h"
#include "root.h"

namespace thames::util::root{

    template<class real>
    real golden_section_search(std::function<real (real)> func, real a, real b, real tol){
        // Declare function for minimisation (root at minimum of absolute of the function)
        std::function<real (real)> f = [func](real x) {return abs(func(x));};

        // Minimise function using the golden section search
        real x0 = thames::util::optimise::golden_section_search(f, a, b, tol);

        // Return root
        return x0;
    }
    template double golden_section_search<double>(std::function<double (double)>, double, double, double);

    template<class real>
    real newton_raphson(const std::function<real (real)> &func, const std::function<real (real)> &dfunc, real xn, real tol){
        // Declare approximation variable
        real xn1;

        // Set converged flag to false
        bool converged = false;

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
    template double newton_raphson<double>(const std::function<double (double)>&, const std::function<double (double)>&, double, double);

}