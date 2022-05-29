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

#include <array>
#include <cmath>
#include <memory>
#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../../include/conversions/dimensional.h"
#include "../../../include/perturbations/geopotential/J2.h"
#include "../../../include/vector/geometry.h"

namespace thames::perturbations::geopotential {

    using thames::conversions::dimensional::DimensionalFactors;

    ///////////
    // Reals //
    ///////////

    template <class T>
    J2<T>::J2(const T& mu, const T& J2, const T& radius, const std::shared_ptr<const DimensionalFactors<T>> factors) : BasePerturbation<T>(factors), m_mu(mu), m_J2(J2), m_radius(radius) {

    }

    template <class T>
    J2<T>::~J2() {

    }

    ////////////
    // Arrays //
    ////////////

    template <class T>
    std::array<T, 3> J2<T>::acceleration_total(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const {
        // Calculate factors
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;
        const T J2 = m_J2;
        const T radius = (m_isNonDimensional) ? m_radius/m_factors->length : m_radius;

        // Extract position components
        const T x = R[0], y = R[1], z = R[2];

        // Calculate range
        const T r = thames::vector::geometry::norm3(R);

        // Precompute factors
        const T J2_fac1 = -1.5*mu*J2*pow(radius, 2.0)/pow(r, 5.0);
        const T J2_fac2 = 5*pow(z, 2.0)/pow(r, 2.0);

        // Declare and calculate perturbing acceleration vector
        const std::array<T, 3> A = {
            J2_fac1*x*(1.0 - J2_fac2),
            J2_fac1*y*(1.0 - J2_fac2),
            J2_fac1*z*(3.0 - J2_fac2)
        };

        // Return perturbing acceleration vector
        return A;
    }

    template <class T>
    T J2<T>::potential(const T& t, const std::array<T, 3>& R) const {
        // Calculate factors
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;
        const T J2 = m_J2;
        const T radius = (m_isNonDimensional) ? m_radius/m_factors->length : m_radius;

        // Extract position components
        const T z = R[2];

        // Calculate range
        const T r = thames::vector::geometry::norm3(R);

        // Calculate cosine of latitude
        const T cphi = z/r;

        // Calculate perturbing potential
        const T U = 0.5*J2*mu/pow(r, 3.0)*pow(radius, 2.0)*(3*pow(cphi, 2.0) - 1.0);

        // Return perturbing potential
        return U;
    }

    /////////////
    // Vectors //
    /////////////

    template <class T>
    std::vector<T> J2<T>::acceleration_total(const T& t, const std::vector<T>& R, const std::vector<T>& V) const {
        // Calculate factors
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;
        const T J2 = m_J2;
        const T radius = (m_isNonDimensional) ? m_radius/m_factors->length : m_radius;

        // Extract position components
        const T x = R[0], y = R[1], z = R[2];

        // Calculate range
        const T r = thames::vector::geometry::norm3(R);

        // Precompute factors
        const T J2_fac1 = -1.5*mu*J2*pow(radius, 2.0)/pow(r, 5.0);
        const T J2_fac2 = 5.0*pow(z, 2.0)/pow(r, 2.0);

        // Declare and calculate perturbing acceleration vector
        const std::vector<T> A = {
            J2_fac1*x*(1.0 - J2_fac2),
            J2_fac1*y*(1.0 - J2_fac2),
            J2_fac1*z*(3.0 - J2_fac2)
        };

        // Return perturbing acceleration vector
        return A;
    }

    template <class T>
    T J2<T>::potential(const T& t, const std::vector<T>& R) const {
        // Calculate factors
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;
        const T J2 = m_J2;
        const T radius = (m_isNonDimensional) ? m_radius/m_factors->length : m_radius;

        // Extract position components
        const T z = R[2];

        // Calculate range
        const T r = thames::vector::geometry::norm3(R);

        // Calculate cosine of latitude
        const T cphi = z/r;

        // Calculate perturbing potential
        const T U = 0.5*J2*mu/pow(r, 3.0)*pow(radius, 2.0)*(3.0*pow(cphi, 2.0) - 1.0);

        // Return perturbing potential
        return U;
    }

    template class J2<double>;

    /////////////////
    // Polynomials //
    /////////////////
    
    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::polynomial;
    using thames::conversions::dimensional::DimensionalFactors;

    template<class T, template<class> class P>
    J2Polynomial<T, P>::J2Polynomial(const T& mu, const T& J2, const T& radius, const std::shared_ptr<const DimensionalFactors<T>> factors) : BasePerturbationPolynomial<T, P>(factors), m_mu(mu), m_J2(J2), m_radius(radius) {

    }

    template<class T, template<class> class P>
    J2Polynomial<T, P>::~J2Polynomial(){

    }

    template<class T, template<class> class P>
    std::vector<P<T>> J2Polynomial<T, P>::acceleration_total(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const {
        // Calculate factors
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;
        const T J2 = m_J2;
        const T radius = (m_isNonDimensional) ? m_radius/m_factors->length : m_radius;

        // Extract position components
        const P<T> x = R[0], y = R[1], z = R[2];

        // Calculate range
        const P<T> r = thames::vector::geometry::norm3(R);

        // Precompute factors
        const P<T> J2_fac1 = -1.5*mu*J2*pow(radius, 2)/pow(r, 5);
        const P<T> J2_fac2 = 5.0*pow(z, 2)/pow(r, 2);

        // Declare and calculate perturbing acceleration vector
        const std::vector<P<T>> A = {
            J2_fac1*x*(1.0 - J2_fac2),
            J2_fac1*y*(1.0 - J2_fac2),
            J2_fac1*z*(3.0 - J2_fac2)
        };

        // Return perturbing acceleration vector
        return A;
    }

    template<class T, template<class> class P>
    P<T> J2Polynomial<T, P>::potential(const T& t, const std::vector<P<T>>& R) const {
        // Calculate factors
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;
        const T J2 = m_J2;
        const T radius = (m_isNonDimensional) ? m_radius/m_factors->length : m_radius;

        // Extract position components
        const P<T> z = R[2];

        // Calculate range
        const P<T> r = thames::vector::geometry::norm3(R);

        // Calculate cosine of latitude
        const P<T> cphi = z/r;

        // Calculate perturbing potential
        const P<T> U = 0.5*J2*mu/pow(r, 3)*pow(radius, 2)*(3.0*pow(cphi, 2) - 1.0);

        // Return perturbing potential
        return U;
    }

    template class J2Polynomial<double, taylor_polynomial>;
    template class J2Polynomial<double, chebyshev_polynomial>;

    #endif

}