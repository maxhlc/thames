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

#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>

#include "../src/thames.h"

int main(int argc, char **argv){
    // Set constants
    double mu = thames::constants::earth::mu;
    double radius = thames::constants::earth::radius;
    double J2 = thames::constants::earth::J2;
    double atol = 1e-13;
    double rtol = 1e-13;

    // Store filepaths as strings
    std::string filepathin(argv[1]), filepathout(argv[2]);

    // Load sample states
    double tstart, tend;
    double tstep = 30;
    int scid;
    thames::constants::statetypes::StateTypes statetype;
    std::vector<std::vector<double>> states;
    thames::io::point::load(filepathin, tstart, tend, scid, statetype, states);

    // Generate polynomials
    std::vector<double> lower, upper;
    std::vector<taylor_polynomial<double>> RVpolynomial, RVpolynomial_propagated;
    thames::conversions::cartesian::cartesian_to_polynomial(states, 4, RVpolynomial, lower, upper);

    // Calculate sample points
    std::vector<std::vector<double>> samples = thames::conversions::cartesian::state_to_sample(states, lower, upper);

    // Non-dimensionalise polynomials
    thames::conversions::dimensional::DimensionalFactors<double> factors;
    thames::conversions::dimensional::cartesian_nondimensionalise(tstart, RVpolynomial, mu, factors);
    radius /= factors.length;
    tend /= factors.time;
    tstep /= factors.time;

    // Declare propagator and perturbation
    thames::perturbations::geopotential::J2Polynomial<double, taylor_polynomial> perturbation(mu, J2, radius);
    thames::propagators::CowellPropagatorPolynomial<double, taylor_polynomial> propagator(mu, &perturbation);

    // Propagate polynomials
    RVpolynomial_propagated = propagator.propagate(tstart, tend, tstep, RVpolynomial, atol, rtol, statetype);

    // Dimensionalise polynomials
    thames::conversions::dimensional::cartesian_dimensionalise(tstart, RVpolynomial_propagated, mu, factors);
    tend *= factors.time;

    // Sample polynomials
    std::vector<std::vector<double>> states_propagated = thames::util::polynomials::evaluate_polynomials(RVpolynomial_propagated, samples);
    
    // Save propagated states
    thames::io::point::save(filepathout, tstart, tend, scid, states_propagated);

    // Return zero
    return 0;
}