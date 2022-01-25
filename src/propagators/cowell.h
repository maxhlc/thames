#ifndef THAMES_PROPAGATORS_COWELL
#define THAMES_PROPAGATORS_COWELL

#include <array>
#include <vector>

#include "../perturbations/baseperturbation.h"

using namespace thames::perturbations::baseperturbation;

namespace thames::propagators::cowell{

    ////////////
    // Arrays //
    ////////////

    /**
     * @brief State derivative for Cowell's method propagation.
     * 
     * @tparam T Numeric type.
     * @param[in] RV Cartesian state.
     * @param[out] RVdot Time derivative of the Cartesian state.
     * @param[in] t Current physical time.
     * @param[in] mu Central body gravitational parameter.
     * @param[in] perturbation Perturbation object.
     */
    template<class T>
    void derivative(const std::array<T, 6>& RV, std::array<T, 6>& RVdot, const T t, const T& mu, const BasePerturbation<T>& perturbation);

    /**
     * @brief Propagate Cartesian state using Cowell's method.
     * 
     * @tparam T Numeric type.
     * @param[in] tstart Propagation start time in physical time.
     * @param[in] tend Propagation end time in physical time.
     * @param[in] tstep Initial timestep for propagation.
     * @param[in] RV Initial Cartesian state.
     * @param[in] mu Central body gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @param[in] atol Solver absolute tolerance.
     * @param[in] rtol Solver relative tolerance.
     * @return std::array<T, 6> Final Cartesian state.
     */
    template<class T>
    std::array<T, 6> propagate(T tstart, T tend, T tstep, std::array<T, 6> RV, T mu, const BasePerturbation<T>& perturbation, T atol = 1e-10, T rtol = 1e-10);

    /////////////
    // Vectors //
    /////////////

    /**
     * @brief State derivative for Cowell's method propagation.
     * 
     * @tparam T Numeric type.
     * @param[in] RV Cartesian state.
     * @param[out] RVdot Time derivative of the Cartesian state.
     * @param[in] t Current physical time.
     * @param[in] mu Central body gravitational parameter.
     * @param[in] perturbation Perturbation object.
     */
    template<class T>
    void derivative(const std::vector<T>& RV, std::vector<T>& RVdot, const T t, const T& mu, const BasePerturbation<T>& perturbation);

    /**
     * @brief Propagate Cartesian state using Cowell's method.
     * 
     * @tparam T Numeric type.
     * @param[in] tstart Propagation start time in physical time.
     * @param[in] tend Propagation end time in physical time.
     * @param[in] tstep Initial timestep for propagation.
     * @param[in] RV Initial Cartesian state.
     * @param[in] mu Central body gravitational parameter.
     * @param[in] perturbation Perturbation object.
     * @param[in] atol Solver absolute tolerance.
     * @param[in] rtol Solver relative tolerance.
     * @return std::vector<T> Final Cartesian state.
     */
    template<class T>
    std::vector<T> propagate(T tstart, T tend, T tstep, std::vector<T> RV, T mu, const BasePerturbation<T>& perturbation, T atol = 1e-10, T rtol = 1e-10);

}

#endif