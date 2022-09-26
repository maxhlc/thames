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

#include <cmath>
#include <stdexcept>

#ifdef THAMES_USE_SMARTUQ
#include "../../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../../include/perturbations/atmosphere/baseatmospheremodel.h"

namespace thames::perturbations::atmosphere::models {

    ///////////
    // Reals //
    ///////////

    template<class T>
    BaseAtmosphereModel<T>::BaseAtmosphereModel() {

    }

    template<class T>
    BaseAtmosphereModel<T>::~BaseAtmosphereModel() {

    }

    template<class T>
    T BaseAtmosphereModel<T>::density(T alt) const {
        // Return zero density
        return 0.0;
    }

    template class BaseAtmosphereModel<double>;

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::polynomial;

    template<class T, template <class> class P>
    BaseAtmosphereModelPolynomial<T, P>::BaseAtmosphereModelPolynomial() {

    }

    template<class T, template <class> class P>
    BaseAtmosphereModelPolynomial<T, P>::~BaseAtmosphereModelPolynomial() {

    }

    template<class T, template <class> class P>
    P<T> BaseAtmosphereModelPolynomial<T, P>::density(P<T> alt) const {
        // Return zero density
        int nvar = alt.get_nvar();
        int degree = alt.get_degree();
        P<T> rho(nvar, degree);
        return rho;
    }

    template class BaseAtmosphereModelPolynomial<double, taylor_polynomial>;
    template class BaseAtmosphereModelPolynomial<double, chebyshev_polynomial>;

    #endif

}