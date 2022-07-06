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

#include <boost/numeric/odeint.hpp>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Integrators/rk4.h"
#include "../../external/smart-uq/include/Integrators/rk45.h"
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../include/constants/statetypes.h"
#include "../../include/conversions/dimensional.h"
#include "../../include/propagators/basepropagator.h"
#include "../../include/propagators/cowell.h"
#include "../../include/perturbations/baseperturbation.h"
#include "../../include/vector/arithmeticoverloads.h"
#include "../../include/vector/geometry.h"

namespace thames::propagators {

    using thames::constants::statetypes::CARTESIAN;
    using thames::perturbations::baseperturbation::BasePerturbation;
    using namespace thames::vector::arithmeticoverloads;
    using thames::conversions::dimensional::DimensionalFactors;

    ///////////
    // Reals //
    ///////////

    template<class T>
    CowellPropagator<T>::CowellPropagator(const T& mu, const std::shared_ptr<BasePerturbation<T>> perturbation, const std::shared_ptr<DimensionalFactors<T>> factors) : BasePropagator<T>(mu, perturbation, factors, CARTESIAN) {

    }

    ////////////
    // Arrays //
    ////////////

    template<class T>
    void CowellPropagator<T>::derivative(const std::array<T, 6>& RV, std::array<T, 6>& RVdot, const T t) const {
        // Calculate factors
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;

        // Extract Cartesian state vectors
        std::array<T, 3> R = {RV[0], RV[1], RV[2]};
        std::array<T, 3> V = {RV[3], RV[4], RV[5]};

        // Calculate range
        T r = thames::vector::geometry::norm3(R);

        // Calculate perturbing acceleration
        std::array<T, 3> F = m_perturbation->acceleration_total(t, R, V);

        // Calculate central body acceleration
        std::array<T, 3> G = -mu/pow(r, 3.0)*R;

        // Calculate acceleration
        std::array<T, 3> A = G + F;

        // Store state derivative
        for(unsigned int ii=0; ii<3; ii++){
            RVdot[ii] = V[ii];
            RVdot[ii+3] = A[ii];
        }
    }

    /////////////
    // Vectors //
    /////////////

    template<class T>
    void CowellPropagator<T>::derivative(const std::vector<T>& RV, std::vector<T>& RVdot, const T t) const {
        // Calculate factors
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;

        // Extract Cartesian state vectors
        std::vector<T> R = {RV[0], RV[1], RV[2]};
        std::vector<T> V = {RV[3], RV[4], RV[5]};

        // Calculate range
        T r = thames::vector::geometry::norm3(R);

        // Calculate perturbing acceleration
        std::vector<T> F = m_perturbation->acceleration_total(t, R, V);

        // Calculate central body acceleration
        std::vector<T> G = -mu/pow(r, 3.0)*R;

        // Calculate acceleration
        std::vector<T> A = G + F;

        // Store state derivative
        for(unsigned int ii=0; ii<3; ii++){
            RVdot[ii] = V[ii];
            RVdot[ii+3] = A[ii];
        }
    }

    template class CowellPropagator<double>;

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::integrator;
    using namespace smartuq::polynomial;
    using thames::perturbations::baseperturbation::BasePerturbationPolynomial;
    using thames::propagators::basepropagator::BasePropagatorPolynomial;
    using thames::propagators::basepropagator::BasePropagatorPolynomialDynamics;

    template<class T, template<class> class P>
    CowellPropagatorPolynomialDynamics<T, P>::CowellPropagatorPolynomialDynamics(const T& mu, const std::shared_ptr<BasePerturbationPolynomial<T, P>> perturbation, const std::shared_ptr<const DimensionalFactors<T>> factors) : BasePropagatorPolynomialDynamics<T, P>("Cowell", mu, perturbation, factors) {

    }

    template<class T, template<class> class P>
    CowellPropagatorPolynomialDynamics<T, P>::~CowellPropagatorPolynomialDynamics() {

    }

    template<class T, template<class> class P>
    int CowellPropagatorPolynomialDynamics<T, P>::evaluate(const T& t, const std::vector<P<T>>& RV, std::vector<P<T>>& RVdot) const {
        // Calculate factors
        const T mu = (m_isNonDimensional) ? m_mu/m_factors->grav : m_mu;

        // Extract Cartesian state vectors
        std::vector<P<T>> R = {RV[0], RV[1], RV[2]};
        std::vector<P<T>> V = {RV[3], RV[4], RV[5]};

        // Calculate range
        P<T> r = thames::vector::geometry::norm3(R);

        // Calculate perturbing acceleration
        std::vector<P<T>> F = m_perturbation->acceleration_total(t, R, V);

        // Calculate central body acceleration
        std::vector<P<T>> G = -mu/pow(r, 3)*R;

        // Calculate acceleration
        std::vector<P<T>> A = G + F;

        // Update derivative
        RVdot = {V[0], V[1], V[2], A[0], A[1], A[2]};

        // Return zero
        return 0;
    }

    template class CowellPropagatorPolynomialDynamics<double, taylor_polynomial>;
    template class CowellPropagatorPolynomialDynamics<double, chebyshev_polynomial>;

    template<class T, template<class> class P>
    CowellPropagatorPolynomial<T, P>::CowellPropagatorPolynomial(const T& mu, const std::shared_ptr<BasePerturbationPolynomial<T, P>> perturbation, const std::shared_ptr<DimensionalFactors<T>> factors) : BasePropagatorPolynomial<T, P>(mu, perturbation, factors, std::make_shared<CowellPropagatorPolynomialDynamics<T, P>>(mu, perturbation, factors), CARTESIAN) {

    }

    template<class T, template<class> class P>
    CowellPropagatorPolynomial<T, P>::~CowellPropagatorPolynomial() {
        
    }

    template class CowellPropagatorPolynomial<double, taylor_polynomial>;
    template class CowellPropagatorPolynomial<double, chebyshev_polynomial>;

    #endif

}