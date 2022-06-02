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

#include "../include/thames.h"

int main(int argc, char **argv){
    // Set constants
    double mu = thames::constants::earth::mu;
    double radius = thames::constants::earth::radius;
    double J2 = thames::constants::earth::J2;
    double atol, rtol;

    // Store filepaths as strings
    std::string filepathin(argv[1]), filepathout(argv[2]);

    // Load sample states
    double tstart, tend, tstep = 30;
    int scid, degree;
    thames::constants::statetypes::StateTypes statetype;
    std::vector<std::vector<double>> states;
    thames::io::point::load(filepathin, tstart, tend, scid, statetype, degree, atol, rtol, states);

    // Create propagator options
    thames::settings::PropagatorParameters<double> options;
    options.absoluteTolerance = atol;
    options.relativeTolerance = rtol;

    // Declare factors
    auto factors = std::make_shared<thames::conversions::dimensional::DimensionalFactors<double>>();

    // Declare perturbation
    auto perturbation = std::make_shared<thames::perturbations::geopotential::J2Polynomial<double, smartuq::polynomial::taylor_polynomial>>(mu, J2, radius, factors);

    // Declare propagators
    thames::propagators::GEqOEPropagatorPolynomial<double, smartuq::polynomial::taylor_polynomial> propagator(mu, perturbation, factors);

    // Propagate states
    std::vector<std::vector<double>> states_propagated = propagator.propagate(tstart, tend, tstep, states, options, statetype, degree);

    // Save propagated states
    thames::io::point::save(filepathout, tstart, tend, scid, statetype, degree, atol, rtol, states_propagated);

    // Return zero
    return 0;
}