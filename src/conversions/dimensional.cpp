/*
MIT License

Copyright (c) 2021-2022 Max Hallgarten La Casta

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <array>
#include <cmath>
#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../include/conversions/dimensional.h"
#include "../../include/conversions/keplerian.h"
#include "../../include/vector/geometry.h"

namespace thames::conversions::dimensional{

    ////////////
    // Arrays //
    ////////////

    template<class T>
    void cartesian_nondimensionalise(T& t, std::array<T, 6>& RV, T& mu, DimensionalFactors<T>& factors){
        // Calculate factors
        factors = calculate_factors(RV, mu);

        // Calculate non-dimensional states
        t /= factors.time;
        mu /= factors.grav;
        RV[0] /= factors.length;
        RV[1] /= factors.length;
        RV[2] /= factors.length;
        RV[3] /= factors.velocity;
        RV[4] /= factors.velocity;
        RV[5] /= factors.velocity;
    }
    template void cartesian_nondimensionalise<double>(double&, std::array<double, 6>&, double&, DimensionalFactors<double>&);

    template<class T>
    void cartesian_dimensionalise(T& t, std::array<T, 6>& RV, T& mu, const DimensionalFactors<T>& factors){
        // Calculate dimensional states
        t *= factors.time;
        mu *= factors.grav;
        RV[0] *= factors.length;
        RV[1] *= factors.length;
        RV[2] *= factors.length;
        RV[3] *= factors.velocity;
        RV[4] *= factors.velocity;
        RV[5] *= factors.velocity;
    }
    template void cartesian_dimensionalise<double>(double&, std::array<double, 6>&, double&, const DimensionalFactors<double>&);

    template<class T>
    DimensionalFactors<T> calculate_factors(const std::array<T, 6> RV, const T mu){
        // Create factors struct
        DimensionalFactors<T> factors;

        // Extract position and velocity vectors
        std::array<T, 3> R = {RV[0], RV[1], RV[2]};
        std::array<T, 3> V = {RV[3], RV[4], RV[5]};

        // Calculate state magnitudes
        T r = thames::vector::geometry::norm3(R);
        T v = thames::vector::geometry::norm3(V);

        // Calculate semi-major axis (length factor)
        factors.length = 1.0/(2.0/r - pow(v, 2.0)/mu);

        // Calculate circular orbit velocity (velocity factor)
        factors.velocity = sqrt(mu/factors.length);

        // Calculate time factor
        factors.time = sqrt(pow(factors.length, 3.0)/mu);

        // Store gravitational parameter factor
        factors.grav = mu;

        // Return factors
        return factors;
    }
    template DimensionalFactors<double> calculate_factors(const std::array<double, 6>, const double);

    /////////////
    // Vectors //
    /////////////

    template<class T>
    void cartesian_nondimensionalise(T& t, std::vector<T>& RV, T& mu, DimensionalFactors<T>& factors){
        // Calculate factors
        factors = calculate_factors(RV, mu);

        // Calculate non-dimensional states
        t /= factors.time;
        mu /= factors.grav;
        RV[0] /= factors.length;
        RV[1] /= factors.length;
        RV[2] /= factors.length;
        RV[3] /= factors.velocity;
        RV[4] /= factors.velocity;
        RV[5] /= factors.velocity;
    }
    template void cartesian_nondimensionalise<double>(double&, std::vector<double>&, double&, DimensionalFactors<double>&);

    template<class T>
    void cartesian_dimensionalise(T &t, std::vector<T>& RV, T& mu, const DimensionalFactors<T>& factors){
        // Calculate dimensional states
        t *= factors.time;
        mu *= factors.grav;
        RV[0] *= factors.length;
        RV[1] *= factors.length;
        RV[2] *= factors.length;
        RV[3] *= factors.velocity;
        RV[4] *= factors.velocity;
        RV[5] *= factors.velocity;
    }
    template void cartesian_dimensionalise<double>(double&, std::vector<double>&, double&, const DimensionalFactors<double>&);

    template<class T>
    DimensionalFactors<T> calculate_factors(const std::vector<T> RV, const T mu){
        // Create factors struct
        DimensionalFactors<T> factors;

        // Extract position and velocity vectors
        std::vector<T> R = {RV[0], RV[1], RV[2]};
        std::vector<T> V = {RV[3], RV[4], RV[5]};

        // Calculate state magnitudes
        T r = thames::vector::geometry::norm3(R);
        T v = thames::vector::geometry::norm3(V);

        // Calculate semi-major axis (length factor)
        factors.length = 1.0/(2.0/r - pow(v, 2.0)/mu);

        // Calculate circular orbit velocity (velocity factor)
        factors.velocity = sqrt(mu/factors.length);

        // Calculate time factor
        factors.time = sqrt(pow(factors.length, 3.0)/mu);

        // Store gravitational parameter factor
        factors.grav = mu;

        // Return factors
        return factors;
    }
    template DimensionalFactors<double> calculate_factors(const std::vector<double>, const double);

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::polynomial;

    template<class T, template<class> class P>
    void cartesian_nondimensionalise(T& t, std::vector<P<T>>& RV, T& mu, DimensionalFactors<T>& factors){
        // Extract central state vector
        std::vector<T> RVcentral = {RV[0].get_coeffs()[0], RV[1].get_coeffs()[0], RV[2].get_coeffs()[0],
                                    RV[0].get_coeffs()[3], RV[4].get_coeffs()[0], RV[5].get_coeffs()[0]};

        // Calculate factors
        factors = calculate_factors(RVcentral, mu);

        // Calculate non-dimensional states
        t /= factors.time;
        mu /= factors.grav;
        RV[0] /= factors.length;
        RV[1] /= factors.length;
        RV[2] /= factors.length;
        RV[3] /= factors.velocity;
        RV[4] /= factors.velocity;
        RV[5] /= factors.velocity;
    }
    template void cartesian_nondimensionalise(double&, std::vector<taylor_polynomial<double>>&, double&, DimensionalFactors<double>&);

    template<class T, template<class> class P>
    void cartesian_dimensionalise(T& t, std::vector<P<T>>& RV, T& mu, DimensionalFactors<T>& factors){
        // Calculate dimensional states
        t *= factors.time;
        mu *= factors.grav;
        RV[0] *= factors.length;
        RV[1] *= factors.length;
        RV[2] *= factors.length;
        RV[3] *= factors.velocity;
        RV[4] *= factors.velocity;
        RV[5] *= factors.velocity;
    }
    template void cartesian_dimensionalise(double&, std::vector<taylor_polynomial<double>>&, double&, DimensionalFactors<double>&);

    #endif

}