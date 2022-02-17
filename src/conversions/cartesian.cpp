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
using namespace smartuq::polynomial;
#endif

#include "cartesian.h"

namespace thames::conversions::cartesian {

    #ifdef THAMES_USE_SMARTUQ

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

    #endif

}