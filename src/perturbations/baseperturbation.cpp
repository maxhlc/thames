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
#include <memory>
#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
#endif

#include "../../include/perturbations/baseperturbation.h"

namespace thames::perturbations::baseperturbation{

    using thames::conversions::dimensional::DimensionalFactors;

    ///////////
    // Reals //
    ///////////

    template<class T>
    BasePerturbation<T>::BasePerturbation(const std::shared_ptr<const DimensionalFactors<T>> factors) : m_factors(factors) {

    };

    template<class T>
    BasePerturbation<T>::~BasePerturbation(){

    };

    template<class T>
    bool BasePerturbation<T>::get_nondimensional() {
        return m_isNonDimensional;
    }

    template<class T>
    void BasePerturbation<T>::set_nondimensional(bool isNonDimensional) {
        m_isNonDimensional = isNonDimensional;
    }

    template<class T>
    std::vector<T> BasePerturbation<T>::acceleration_total(const T& t, const std::vector<T>& R, const std::vector<T>& V) const{
        std::vector<T> F = {0.0, 0.0, 0.0};
        return F;
    };

    template<class T>
    std::vector<T> BasePerturbation<T>::acceleration_nonpotential(const T& t, const std::vector<T>& R, const std::vector<T>& V) const{
        std::vector<T> F = {0.0, 0.0, 0.0};
        return F;
    };

    template<class T>
    T BasePerturbation<T>::potential(const T& t, const std::vector<T>& R) const{
        T U = 0.0;
        return U;
    }

    template<class T>
    T BasePerturbation<T>::potential_derivative(const T& t, const std::vector<T>& R, const std::vector<T>& V) const{
        T Ut = 0.0;
        return Ut;
    }

    template class BasePerturbation<double>;

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    using namespace smartuq::polynomial;
    using thames::conversions::dimensional::DimensionalFactors;
    
    template<class T, template<class> class P>
    BasePerturbationPolynomial<T, P>::BasePerturbationPolynomial(const std::shared_ptr<const DimensionalFactors<T>> factors) : m_factors(factors) {

    };

    template<class T, template<class> class P>
    BasePerturbationPolynomial<T, P>::~BasePerturbationPolynomial(){

    };

    template<class T, template<class> class P>
    bool BasePerturbationPolynomial<T, P>::get_nondimensional() {
        return m_isNonDimensional;
    }

    template<class T, template<class> class P>
    void BasePerturbationPolynomial<T, P>::set_nondimensional(const bool isNonDimensional) {
        m_isNonDimensional = isNonDimensional;
    }

    template<class T, template<class> class P>
    std::vector<P<T>> BasePerturbationPolynomial<T, P>::acceleration_total(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const {
        int nvar = R[0].get_nvar();
        int degree = R[0].get_degree();
        P<T> poly(nvar, degree);
        std::vector<P<T>> F = {poly, poly, poly};
        return F;
    };

    template<class T, template<class> class P>
    std::vector<P<T>> BasePerturbationPolynomial<T, P>::acceleration_nonpotential(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const {
        int nvar = R[0].get_nvar();
        int degree = R[0].get_degree();
        P<T> poly(nvar, degree);
        std::vector<P<T>> F = {poly, poly, poly};
        return F;
    };

    template<class T, template<class> class P>
    P<T> BasePerturbationPolynomial<T, P>::potential(const T& t, const std::vector<P<T>>& R) const {
        int nvar = R[0].get_nvar();
        int degree = R[0].get_degree();
        P<T> U(nvar, degree);
        return U;
    }

    template<class T, template<class> class P>
    P<T> BasePerturbationPolynomial<T, P>::potential_derivative(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const {
        int nvar = R[0].get_nvar();
        int degree = R[0].get_degree();
        P<T> Ut(nvar, degree);
        return Ut;
    }

    template class BasePerturbationPolynomial<double, taylor_polynomial>;
    template class BasePerturbationPolynomial<double, chebyshev_polynomial>;

    #endif

}