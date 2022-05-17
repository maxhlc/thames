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

#include "../../../include/perturbations/atmosphere/atmospheremodel.h"
#include "../../../include/perturbations/atmosphere/drag.h"
#include "../../../include/perturbations/atmosphere/ussa76.h"
#include "../../../include/perturbations/baseperturbation.h"
#include "../../../include/vector/arithmeticoverloads.h"
#include "../../../include/vector/geometry.h"

namespace thames::perturbations::atmosphere::drag {

    using thames::perturbations::atmosphere::models::AtmosphereModel;
    using thames::perturbations::baseperturbation::BasePerturbation;

    using namespace thames::vector::arithmeticoverloads;

    template<class T>
    Drag<T>::Drag(const T& radius, const T& w, const T& Cd, const T& A, const T& m, const AtmosphereModels& model, const DimensionalFactors<T>* factors) : BasePerturbation<T>(factors), m_radius(radius), m_w(w), m_Cd(Cd), m_A(A), m_m(m), m_model(AtmosphereModel<T>(model)) {

    }

    template<class T>
    Drag<T>::~Drag() {

    }

    template<class T>
    std::array<T, 3> Drag<T>::acceleration_total(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const {
        return acceleration_nonpotential(t, R, V);
    }

    template<class T>
    std::array<T, 3> Drag<T>::acceleration_nonpotential(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const {
        // Calculate factors
        T radius = (m_isNonDimensional) ? m_radius/(m_factors->length) : m_radius;
        T w = (m_isNonDimensional) ? m_w*(m_factors->time) : m_w;
        T Cd = m_Cd;
        T A = (m_isNonDimensional) ? m_A/std::pow(m_factors->length, 2) : m_A;

        // Calculate altitude
        T r = thames::vector::geometry::norm3(R);
        T alt = r - radius;

        // Calculate atmospheric density
        if (m_isNonDimensional)
            alt *= m_factors->length;
        T rho = m_model.density(alt);
        
        // Calculate factors which include mass (cancelled via rho/mass) and non-dimensionalise as required
        T massfac = rho/m_m;
        if (m_isNonDimensional)
            massfac *= std::pow(m_factors->length, 3);

        // Calculate velocity relative to the atmosphere
        std::array<T, 3> W = {0, 0, w};
        std::array<T, 3> Vrel = V - thames::vector::geometry::cross3(W, R);
        T vrel = thames::vector::geometry::norm3(Vrel);
        std::array<T, 3> uv = Vrel/vrel;

        // Calculate acceleration due to drag
        std::array<T, 3> Ad = -0.5*Cd*A*massfac*pow(vrel, 2)*uv;

        // Return acceleration
        return Ad;
    }

    template<class T>
    std::vector<T> Drag<T>::acceleration_total(const T& t, const std::vector<T>& R, const std::vector<T>& V) const {
        return acceleration_nonpotential(t, R, V);
    }

    template<class T>
    std::vector<T> Drag<T>::acceleration_nonpotential(const T& t, const std::vector<T>& R, const std::vector<T>& V) const {
        // Calculate factors
        T radius = (m_isNonDimensional) ? m_radius/(m_factors->length) : m_radius;
        T w = (m_isNonDimensional) ? m_w*(m_factors->time) : m_w;
        T Cd = m_Cd;
        T A = (m_isNonDimensional) ? m_A/std::pow(m_factors->length, 2) : m_A;

        // Calculate altitude
        T r = thames::vector::geometry::norm3(R);
        T alt = r - radius;

        // Calculate atmospheric density
        if (m_isNonDimensional)
            alt *= m_factors->length;
        T rho = m_model.density(alt);
        
        // Calculate factors which include mass (cancelled via rho/mass) and non-dimensionalise as required
        T massfac = rho/m_m;
        if (m_isNonDimensional)
            massfac *= std::pow(m_factors->length, 3);

        // Calculate velocity relative to the atmosphere
        std::vector<T> W = {0, 0, w};
        std::vector<T> Vrel = V - thames::vector::geometry::cross3(W, R);
        T vrel = thames::vector::geometry::norm3(Vrel);
        std::vector<T> uv = Vrel/vrel;

        // Calculate acceleration due to drag
        std::vector<T> Ad = -0.5*Cd*A*massfac*pow(vrel, 2)*uv;

        // Return acceleration
        return Ad;
    }

    template class Drag<double>;

}