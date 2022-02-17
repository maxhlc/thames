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

#include "polynomials.h"

namespace thames::util::polynomials {

    #ifdef THAMES_USE_SMARTUQ

    template<class T, template<class> class P>
    std::vector<T> evaluate_polynomials(const std::vector<P<T>>& polynomials, const std::vector<T>& x) {
        // Declare point state vector
        std::vector<T> numeric;

        // Evaluate each polynomial
        for(std::size_t ii=0; ii<polynomials.size(); ii++)
            numeric.push_back(polynomials[ii].evaluate(x));

        // Return point state vector
        return numeric;
    }
    template std::vector<double> evaluate_polynomials(const std::vector<taylor_polynomial<double>>&, const std::vector<double>&);
    template std::vector<double> evaluate_polynomials(const std::vector<chebyshev_polynomial<double>>&, const std::vector<double>&);

    #endif

}