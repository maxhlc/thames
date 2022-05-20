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

#include "../../include/perturbations/baseperturbation.h"
#include "../../include/perturbations/perturbationcombiner.h"
#include "../../include/vector/arithmeticoverloads.h"

namespace thames::perturbations::perturbationcombiner {

    using thames::conversions::dimensional::DimensionalFactors;
    using thames::perturbations::baseperturbation::BasePerturbation;

    using namespace thames::vector::arithmeticoverloads;

    template<class T>
    PerturbationCombiner<T>::PerturbationCombiner(const DimensionalFactors<T>* const factors) : BasePerturbation<T>(factors) {
        // Ensure all underlying models have same non-dimensional flag
        set_nondimensional(m_isNonDimensional);
    }

    template<class T>
    PerturbationCombiner<T>::PerturbationCombiner(const std::vector<std::shared_ptr<BasePerturbation<T>>>& models, const DimensionalFactors<T>* const factors) : BasePerturbation<T>(factors), m_models(models) {
        // Ensure all underlying models have same non-dimensional flag
        set_nondimensional(m_isNonDimensional);
    }

    template<class T>
    PerturbationCombiner<T>::~PerturbationCombiner() {

    }

    template<class T>
    void PerturbationCombiner<T>::set_nondimensional(const bool isNonDimensional) {
        // Iterate through underlying models to set non-dimensional flag
        for (std::shared_ptr<BasePerturbation<T>> model : m_models)
            model->set_nondimensional(isNonDimensional);
    }

    template<class T>
    void PerturbationCombiner<T>::add_model(const std::shared_ptr<BasePerturbation<T>>& model) {
        // Ensure model has same non-dimensional flag
        model->set_nondimensional(m_isNonDimensional);

        // Add model to class vector
        m_models.push_back(model);
    }

    ////////////
    // Arrays //
    ////////////

    template<class T>
    std::array<T, 3> PerturbationCombiner<T>::acceleration_total(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const {
        /// Declare zero total acceleration
        std::array<T, 3> F = {0.0, 0.0, 0.0};

        // Iterate through underlying models to add to the total acceleration
        for (auto model : m_models)
            F = F + model->acceleration_total(t, R, V);

        // Return acceleration
        return F;
    }

    template<class T>
    std::array<T, 3> PerturbationCombiner<T>::acceleration_nonpotential(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const {
        /// Declare zero non-potential acceleration
        std::array<T, 3> F = {0.0, 0.0, 0.0};

        // Iterate through underlying models to add to the non-potential acceleration
        for (auto model : m_models)
            F = F + model->acceleration_nonpotential(t, R, V);

        // Return acceleration
        return F;
    }

    template<class T>
    T PerturbationCombiner<T>::potential(const T& t, const std::array<T, 3>& R) const {
        /// Declare zero potential
        T U = 0.0;

        // Iterate through underlying models to add to the potential
        for (auto model : m_models)
            U += model->potential(t, R);

        // Return potential
        return U;
    }

    template<class T>
    T PerturbationCombiner<T>::potential_derivative(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const {
        /// Declare zero potential derivative
        T Ut = 0.0;

        // Iterate through underlying models to add to the potential derivative
        for (auto model : m_models)
            Ut += model->potential_derivative(t, R, V);

        // Return potential derivative
        return Ut;
    }

    /////////////
    // Vectors //
    /////////////

    template<class T>
    std::vector<T> PerturbationCombiner<T>::acceleration_total(const T& t, const std::vector<T>& R, const std::vector<T>& V) const {
        /// Declare zero total acceleration
        std::vector<T> F = {0.0, 0.0, 0.0};

        // Iterate through underlying models to add to the total acceleration
        for (auto model : m_models)
            F = F + model->acceleration_total(t, R, V);

        // Return acceleration
        return F;
    }

    template<class T>
    std::vector<T> PerturbationCombiner<T>::acceleration_nonpotential(const T& t, const std::vector<T>& R, const std::vector<T>& V) const {
        /// Declare zero non-potential acceleration
        std::vector<T> F = {0.0, 0.0, 0.0};

        // Iterate through underlying models to add to the non-potential acceleration
        for (auto model : m_models)
            F = F + model->acceleration_nonpotential(t, R, V);

        // Return acceleration
        return F;
    }

    template<class T>
    T PerturbationCombiner<T>::potential(const T& t, const std::vector<T>& R) const {
        /// Declare zero potential
        T U = 0.0;

        // Iterate through underlying models to add to the potential
        for (auto model : m_models)
            U += model->potential(t, R);

        // Return potential
        return U;
    }

    template<class T>
    T PerturbationCombiner<T>::potential_derivative(const T& t, const std::vector<T>& R, const std::vector<T>& V) const {
        /// Declare zero potential derivative
        T Ut = 0.0;

        // Iterate through underlying models to add to the potential derivative
        for (auto model : m_models)
            Ut += model->potential_derivative(t, R, V);

        // Return potential derivative
        return Ut;
    }

    template class PerturbationCombiner<double>;
    
}