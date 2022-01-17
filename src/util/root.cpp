#include <cmath>
#include <functional>

#include "optimise.h"
#include "root.h"

namespace thames::util::root{

    template<class T>
    T golden_section_search(std::function<T (T)> func, T a, T b, T tol){
        // Declare function for minimisation (root at minimum of absolute of the function)
        std::function<T (T)> f = [func](T x) {return fabs(func(x));};

        // Minimise function using the golden section search
        T x0 = thames::util::optimise::golden_section_search(f, a, b, tol);

        // Return root
        return x0;
    }
    template double golden_section_search<double>(std::function<double (double)>, double, double, double);

    template<class T>
    T newton_raphson(const std::function<T (T)> &func, const std::function<T (T)> &dfunc, T xn, T tol){
        // Declare approximation variable
        T xn1;

        // Set converged flag to false
        bool converged = false;

        // Iterate until converged
        while(!converged){
            // Update approximation
            xn1 = xn - func(xn)/dfunc(xn);

            // Converged if update is smaller than tolerance
            if(fabs(xn1 - xn) < tol)
                converged = true;

            // Update previous approximation
            xn = xn1;
        }

        // Return root
        return xn1;
    }
    template double newton_raphson<double>(const std::function<double (double)>&, const std::function<double (double)>&, double, double);

}