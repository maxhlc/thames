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

#include "../../../include/perturbations/atmosphere/wertzp.h"

namespace thames::perturbations::atmosphere::models {

    ///////////
    // Reals //
    ///////////

    template<class T>
    WertzPAtmosphereModel<T>::WertzPAtmosphereModel() {

    }

    template<class T>
    WertzPAtmosphereModel<T>::~WertzPAtmosphereModel() {

    }

    template<class T>
    T WertzPAtmosphereModel<T>::density(T alt) const {
        // Scale altitude
        const T altscaled = 2.0*(alt - m_domain[0])/(m_domain[1] - m_domain[0]) - 1.0;

        // Evaluate polynomial
        T rho = 0.0;
        for (std::size_t ii=0; ii<m_coeff.size(); ii++) {
            rho += m_coeff[ii]*pow(altscaled, ii);
        }
        rho = exp(rho);

        // Return density
        return rho;
    }

    template class WertzPAtmosphereModel<double>;

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::polynomial;

    template<class T, template <class> class P>
    WertzPAtmosphereModelPolynomial<T, P>::WertzPAtmosphereModelPolynomial() {

    }

    template<class T, template <class> class P>
    WertzPAtmosphereModelPolynomial<T, P>::~WertzPAtmosphereModelPolynomial() {

    }

    
    template<class T, template <class> class P>
    P<T> WertzPAtmosphereModelPolynomial<T, P>::density(P<T> alt) const {
        // Scale altitude
        const P<T> altscaled = 2.0*(alt - m_domain[0])/(m_domain[1] - m_domain[0]) - 1.0;

        // Evaluate polynomial
        P<T> rho(alt.get_nvar(), alt.get_degree());
        for (std::size_t ii=0; ii<m_coeff.size(); ii++) {
            rho += m_coeff[ii]*pow(altscaled, ii);
        }
        rho = exp(rho);

        // Return density
        return rho;
    }

    template class WertzPAtmosphereModelPolynomial<double, taylor_polynomial>;
    template class WertzPAtmosphereModelPolynomial<double, chebyshev_polynomial>;

    #endif

}