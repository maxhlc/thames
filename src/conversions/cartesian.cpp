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

#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../include/conversions/cartesian.h"

namespace thames::conversions::cartesian {

    template<class T>
    void states_to_bounds(const std::vector<std::vector<T>>& states, std::vector<T>& lower, std::vector<T>& upper){
        // Resize vectors for minimum and maximum values
        lower.resize(6);
        upper.resize(6);

        // Iterate through states
        for(std::size_t ii=0; ii<states.size(); ii++){
            // Store values for first state, otherwise compare with current minimum/maximum
            if(ii == 0){
                lower = states[ii];
                upper = states[ii];
            } else {
                // Iterate through state variables
                for(unsigned int jj = 0; jj<6; jj++){
                    // Update minimum if state variable is smaller
                    if(lower[jj] > states[ii][jj])
                        lower[jj] = states[ii][jj];

                    // Update maximum if state variable is larger
                    if(upper[jj] < states[ii][jj])
                        upper[jj] = states[ii][jj];                
                }
            }
        }
    }
    template void states_to_bounds(const std::vector<std::vector<double>>& states, std::vector<double>& lower, std::vector<double>& upper);

    template<class T>
    std::vector<T> state_to_sample(const std::vector<T>& state, const std::vector<T>& lower, const std::vector<T>& upper){
        // Declare vector for sample
        std::vector<T> sample(state.size());

        // Iterate through state
        for(std::size_t ii=0; ii<state.size(); ii++){
            // Set sample value to zero when the state has zero range, or scale state value between [-1,1] otherwise
            if((upper[ii] - lower[ii]) == 0.0){
                sample[ii] = 0.0;
            } else {
                sample[ii] = 2.0*(state[ii] - lower[ii])/(upper[ii] - lower[ii]) - 1.0;
            }
        }

        // Return sample point
        return sample;
    }
    template std::vector<double> state_to_sample(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&);

    template<class T>
    std::vector<std::vector<T>> state_to_sample(const std::vector<std::vector<T>>& states, const std::vector<T>& lower, const std::vector<T>& upper){
        // Declare vector for samples
        std::vector<std::vector<T>> samples(states.size(), std::vector<T>(states[0].size()));

        // Iterate through states
        for(std::size_t ii=0; ii<states.size(); ii++)
            samples[ii] = state_to_sample(states[ii], lower, upper);

        // Return samples
        return samples;
    }
    template std::vector<std::vector<double>> state_to_sample(const std::vector<std::vector<double>>&, const std::vector<double>&, const std::vector<double>&);

    #ifdef THAMES_USE_SMARTUQ
    
    using namespace smartuq::polynomial;

    template<class T, template<class> class P>
    void cartesian_to_polynomial(const std::vector<T>& RV, const std::vector<T>& RVunc, int degree, std::vector<P<T>>& RVPolynomial) {
        // Clear polynomial vector
        RVPolynomial.clear();

        // Iterate through Cartesian state variables
        for(unsigned int ii=0; ii<6; ii++){
            // Generate polynomial
            RVPolynomial.push_back(P<T>(6, degree, ii, RV[ii]-RVunc[ii], RV[ii]+RVunc[ii]));

            // Initialise fast multiplication
            RVPolynomial[ii].initialize_M(6, degree);  
        }
    }
    template void cartesian_to_polynomial(const std::vector<double>& RV, const std::vector<double>& RVunc, int degree, std::vector<taylor_polynomial<double>>& RVPolynomial);
    template void cartesian_to_polynomial(const std::vector<double>& RV, const std::vector<double>& RVunc, int degree, std::vector<chebyshev_polynomial<double>>& RVPolynomial);

    template<class T, template<class> class P>
    void cartesian_to_polynomial(const std::vector<std::vector<T>>& RVs, int degree, std::vector<P<T>>& RVPolynomial, std::vector<T>& lower, std::vector<T>& upper) {
        // Clear polynomial vector
        RVPolynomial.clear();

        // Calculate bounds of Cartesian states
        states_to_bounds(RVs, lower, upper);

        // Iterate through Cartesian state variables
        for(unsigned int ii=0; ii<6; ii++){
            // Generate polynomial
            RVPolynomial.push_back(P<T>(6, degree, ii, lower[ii], upper[ii]));

            // Initialise fast multiplication
            RVPolynomial[ii].initialize_M(6, degree);  
        }
    }
    template void cartesian_to_polynomial(const std::vector<std::vector<double>>&, int, std::vector<taylor_polynomial<double>>&, std::vector<double>&, std::vector<double>&);
    template void cartesian_to_polynomial(const std::vector<std::vector<double>>&, int, std::vector<chebyshev_polynomial<double>>&, std::vector<double>&, std::vector<double>&);

    #endif

}