#include <array>
#include <iostream>
#include <iomanip>
#include <vector>

#include "constants/constants.h"
#include "conversions/conversions.h"
#include "perturbations/perturbations.h"
#include "propagators/propagators.h"

int main(){
    std::array<double, 6> RV;
    RV = {6.916000000000002E+03, 0.000000000000000E+00, 0.000000000000000E+00, 0.000000000000000E+00, 4.483946567026534E+00, 8.954234385325996E+00};
    double mu = 3.986004414498200E+05; //thames::constants::earth::mu;
    double tstart = 0.0;
    double tstep = 60.0;
    double tend = 86400.0 * 365.25;
    double radius = 6.378136460000000E+03; //thames::constants::earth::radius;
    double J2 = 1.082626111e-3; //thames::constants::earth::J2;

    thames::perturbations::geopotential::J2<double> perturbation(mu, J2, radius);

    std::array<double, 6> state_prop = thames::propagators::cowell::propagate<double>(tstart, tend, tstep, RV, mu, perturbation, 1e-13, 1e-13);

    std::cout << std::setprecision(16);
    for(unsigned int ii=0; ii<6; ii++)
        std::cout << state_prop[ii] << "\n";

    std::cout << "THAMES executed successfully" << std::endl;

    return 0;
}