#include <array>
#include <cmath>
#include <vector>

#include "J2.h"
#include "../../vector/geometry.h"

namespace thames::perturbations::geopotential{

    template <class T>
    J2<T>::J2(T mu, T J2, T radius){
        m_mu = mu;
        m_J2 = J2;
        m_radius = radius;
    }

    template <class T>
    std::array<T, 3> J2<T>::acceleration_total(T t, std::array<T, 3> R, std::array<T, 3> V) const{
        // Extract position components
        T x = R[0], y = R[1], z = R[2];

        // Calculate range
        T r = thames::vector::geometry::norm3<T>(R);

        // Precompute factors
        T J2_fac1 = -1.5*m_mu*m_J2*pow(m_radius, 2.0)/pow(r, 5.0);
        T J2_fac2 = 5.0*pow(z, 2.0)/pow(r, 2.0);

        // Declare and calculate perturbing acceleration vector
        std::array<T, 3> A;
        A[0] =  J2_fac1*x*(1.0 - J2_fac2),
        A[1] =  J2_fac1*y*(1.0 - J2_fac2),
        A[2] =  J2_fac1*z*(3.0 - J2_fac2);

        // Return perturbing acceleration vector
        return A;
    }

    template <class T>
    std::vector<T> J2<T>::acceleration_total(T t, std::vector<T> R, std::vector<T> V) const{
        // Extract position components
        T x = R[0], y = R[1], z = R[2];

        // Calculate range
        T r = thames::vector::geometry::norm3<T>(R);

        // Precompute factors
        T J2_fac1 = -1.5*m_mu*m_J2*pow(m_radius, 2.0)/pow(r, 5.0);
        T J2_fac2 = 5.0*pow(z, 2.0)/pow(r, 2.0);

        // Declare and calculate perturbing acceleration vector
        std::vector<T> A(3);
        A[0] =  J2_fac1*x*(1.0 - J2_fac2),
        A[1] =  J2_fac1*y*(1.0 - J2_fac2),
        A[2] =  J2_fac1*z*(3.0 - J2_fac2);

        // Return perturbing acceleration vector
        return A;
    }

    template <class T>
    T J2<T>::potential(T t, std::array<T, 3> R) const{
        // Extract position components
        T z = R[2];

        // Calculate range
        T r = thames::vector::geometry::norm3<T>(R);

        // Calculate cosine of latitude
        T cphi = z/r;

        // Calculate perturbing potential
        T U = 0.5*m_J2*m_mu/pow(r, 3.0)*pow(m_radius, 2.0)*(3.0*pow(cphi, 2.0) - 1.0);

        // Return perturbing potential
        return U;
    }

    template <class T>
    T J2<T>::potential(T t, std::vector<T> R) const{
        // Extract position components
        T z = R[2];

        // Calculate range
        T r = thames::vector::geometry::norm3<T>(R);

        // Calculate cosine of latitude
        T cphi = z/r;

        // Calculate perturbing potential
        T U = 0.5*m_J2*m_mu/pow(r, 3.0)*pow(m_radius, 2.0)*(3.0*pow(cphi, 2.0) - 1.0);

        // Return perturbing potential
        return U;
    }

    template class J2<double>;

}