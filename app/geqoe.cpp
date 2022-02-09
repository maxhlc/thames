#include <array>
#include <iostream>
#include <iomanip>
#include <vector>

#include "../src/thames.h"

int main(){
    // Set parameters
    std::array<double, 6> RV;
    RV = {6.916000000000002E+03, 0.000000000000000E+00, 0.000000000000000E+00, 0.000000000000000E+00, 4.483946567026534E+00, 8.954234385325996E+00};
    double mu = 3.986004414498200E+05;
    double tstart = 0.0;
    double tstep = 60.0;
    double tend = 86400.0 * 365.25;
    double radius = 6.378136460000000E+03;
    double J2 = 1.082626111e-3;

    // Non-dimensionalise state
    thames::conversions::dimensional::DimensionalFactors<double> factors;
    thames::conversions::dimensional::cartesian_nondimensionalise(tstart, RV, mu, factors);
    tstep /= factors.time;
    tend /= factors.time;
    radius /= factors.length;

    // Declare J2 perturbation
    thames::perturbations::geopotential::J2<double> perturbation(mu, J2, radius);

    // Declare propagator
    thames::propagators::GEqOEPropagator<double> propagator(mu, &perturbation);

    // Propagate state
    std::array<double, 6> state_prop = propagator.propagate(tstart, tend, tstep, RV, 1e-13, 1e-13);

    // Dimensionalise state
    thames::conversions::dimensional::cartesian_dimensionalise(tend, state_prop, mu, factors);

    // Output final state
    std::cout << std::setprecision(16);
    std::cout << "GEqOE point propagation:\n";
    for(unsigned int ii=0; ii<6; ii++)
        std::cout << state_prop[ii] << "\n";
    std::cout << "THAMES executed successfully" << std::endl;

    return 0;
}