#include <array>
#include <cmath>
#include <vector>

#ifdef THAMES_USE_SMARTUQ
#include "../../../external/smart-uq/include/Polynomial/smartuq_polynomial.h"
using namespace smartuq::polynomial;
#endif

#include "J2.h"
#include "../../vector/geometry.h"

namespace thames::perturbations::geopotential{

    ///////////
    // Reals //
    ///////////

    template <class T>
    J2<T>::J2(const T& mu, const T& J2, const T& radius) : m_mu(mu), m_J2(J2), m_radius(radius) {

    }

    template <class T>
    J2<T>::~J2() {

    }

    ////////////
    // Arrays //
    ////////////

    template <class T>
    std::array<T, 3> J2<T>::acceleration_total(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const{
        // Extract position components
        T x = R[0], y = R[1], z = R[2];

        // Calculate range
        T r = thames::vector::geometry::norm3(R);

        // Precompute factors
        T J2_fac1 = -1.5*m_mu*m_J2*pow(m_radius, 2.0)/pow(r, 5.0);
        T J2_fac2 = 5.0*pow(z, 2.0)/pow(r, 2.0);

        // Declare and calculate perturbing acceleration vector
        std::array<T, 3> A = {
            J2_fac1*x*(1.0 - J2_fac2),
            J2_fac1*y*(1.0 - J2_fac2),
            J2_fac1*z*(3.0 - J2_fac2)
        };

        // Return perturbing acceleration vector
        return A;
    }

    template <class T>
    T J2<T>::potential(const T& t, const std::array<T, 3>& R) const{
        // Extract position components
        T z = R[2];

        // Calculate range
        T r = thames::vector::geometry::norm3(R);

        // Calculate cosine of latitude
        T cphi = z/r;

        // Calculate perturbing potential
        T U = 0.5*m_J2*m_mu/pow(r, 3.0)*pow(m_radius, 2.0)*(3.0*pow(cphi, 2.0) - 1.0);

        // Return perturbing potential
        return U;
    }

    /////////////
    // Vectors //
    /////////////

    template <class T>
    std::vector<T> J2<T>::acceleration_total(const T& t, const std::vector<T>& R, const std::vector<T>& V) const{
        // Extract position components
        T x = R[0], y = R[1], z = R[2];

        // Calculate range
        T r = thames::vector::geometry::norm3(R);

        // Precompute factors
        T J2_fac1 = -1.5*m_mu*m_J2*pow(m_radius, 2.0)/pow(r, 5.0);
        T J2_fac2 = 5.0*pow(z, 2.0)/pow(r, 2.0);

        // Declare and calculate perturbing acceleration vector
        std::vector<T> A = {
            J2_fac1*x*(1.0 - J2_fac2),
            J2_fac1*y*(1.0 - J2_fac2),
            J2_fac1*z*(3.0 - J2_fac2)
        };

        // Return perturbing acceleration vector
        return A;
    }

    template <class T>
    T J2<T>::potential(const T& t, const std::vector<T>& R) const{
        // Extract position components
        T z = R[2];

        // Calculate range
        T r = thames::vector::geometry::norm3(R);

        // Calculate cosine of latitude
        T cphi = z/r;

        // Calculate perturbing potential
        T U = 0.5*m_J2*m_mu/pow(r, 3.0)*pow(m_radius, 2.0)*(3.0*pow(cphi, 2.0) - 1.0);

        // Return perturbing potential
        return U;
    }

    template class J2<double>;

    /////////////////
    // Polynomials //
    /////////////////
    
    #ifdef THAMES_USE_SMARTUQ

    template<class T, template<class> class P>
    J2Polynomial<T, P>::J2Polynomial(const T& mu, const T& J2, const T& radius) : m_mu(mu), m_J2(J2), m_radius(radius) {

    }

    template<class T, template<class> class P>
    J2Polynomial<T, P>::~J2Polynomial(){

    }

    template<class T, template<class> class P>
    std::vector<P<T>> J2Polynomial<T, P>::acceleration_total(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const {
        // Extract position components
        P<T> x = R[0], y = R[1], z = R[2];

        // Calculate range
        P<T> r = thames::vector::geometry::norm3(R);

        // Precompute factors
        P<T> J2_fac1 = -1.5*m_mu*m_J2*pow(m_radius, 2)/pow(r, 5);
        P<T> J2_fac2 = 5.0*pow(z, 2)/pow(r, 2);

        // Declare and calculate perturbing acceleration vector
        std::vector<P<T>> A = {
            J2_fac1*x*(1.0 - J2_fac2),
            J2_fac1*y*(1.0 - J2_fac2),
            J2_fac1*z*(3.0 - J2_fac2)
        };

        // Return perturbing acceleration vector
        return A;
    }

    template<class T, template<class> class P>
    P<T> J2Polynomial<T, P>::potential(const T& t, const std::vector<P<T>>& R) const {
        // Extract position components
        P<T> z = R[2];

        // Calculate range
        P<T> r = thames::vector::geometry::norm3(R);

        // Calculate cosine of latitude
        P<T> cphi = z/r;

        // Calculate perturbing potential
        P<T> U = 0.5*m_J2*m_mu/pow(r, 3)*pow(m_radius, 2)*(3.0*pow(cphi, 2) - 1.0);

        // Return perturbing potential
        return U;
    }

    template class J2Polynomial<double, taylor_polynomial>;
    template class J2Polynomial<double, chebyshev_polynomial>;

    #endif

}