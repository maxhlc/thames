#include <array>
#include <cmath>
#include <vector>

#include "dimensional.h"
#include "keplerian.h"
#include "../vector/geometry.h"

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

}