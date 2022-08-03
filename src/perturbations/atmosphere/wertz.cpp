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

#include "../../../include/perturbations/atmosphere/wertz.h"

namespace thames::perturbations::atmosphere::models {

    ///////////
    // Reals //
    ///////////

    template<class T>
    WertzAtmosphereModel<T>::WertzAtmosphereModel() {

    }

    template<class T>
    WertzAtmosphereModel<T>::~WertzAtmosphereModel() {

    }

    template<class T>
    T WertzAtmosphereModel<T>::density(T alt) const {
        // Declare index variable
        std::size_t ii;

        // Handle altitudes outside of the range
        T altselect = alt;
        if (altselect > m_geo.back()) {
            altselect = m_geo.back();
        } else if (altselect < m_geo.front()) {
            altselect = m_geo.front();
        }

        // Determine interpolation interval
        for (std::size_t jj = 0; jj < m_geo.size() - 1; jj++) {
            if (altselect >= m_geo[jj] && altselect < m_geo[jj+1])
                ii = jj;
        }
        if (altselect >= m_geo.back())
            ii = m_geo.size() - 1;

        // Exponential interpolation
        T rho = m_rho[ii]*std::exp(-(alt - m_geo[ii])/m_scale[ii]);

        // Return density
        return rho;
    }

    template class WertzAtmosphereModel<double>;

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::polynomial;

    template<class T, template <class> class P>
    WertzAtmosphereModelPolynomial<T, P>::WertzAtmosphereModelPolynomial() {

    }

    template<class T, template <class> class P>
    WertzAtmosphereModelPolynomial<T, P>::~WertzAtmosphereModelPolynomial() {

    }

    
    template<class T, template <class> class P>
    P<T> WertzAtmosphereModelPolynomial<T, P>::density(P<T> alt) const {
        throw std::runtime_error("Density profile not implemented");
    }

    template class WertzAtmosphereModelPolynomial<double, taylor_polynomial>;
    template class WertzAtmosphereModelPolynomial<double, chebyshev_polynomial>;

    #endif

}