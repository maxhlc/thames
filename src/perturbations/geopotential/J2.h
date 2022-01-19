#ifndef THAMES_PERTURBATIONS_GEOPOTENTIAL_J2
#define THAMES_PERTURBATIONS_GEOPOTENTIAL_J2

#include <array>
#include <vector>

#include "../baseperturbation.h"

using namespace thames::perturbations::baseperturbation;

namespace thames::perturbations::geopotential{

    /**
     * @brief Class for the perturbation resulting from the J2-term.
     * 
     * @tparam T Numeric type.
     */
    template<class T>
    class J2 : public BasePerturbation<T> {
        private:

            /// Central body gravitational parameter.
            T m_mu;

            /// Central body J2-term.
            T m_J2;

            /// Central body radius.
            T m_radius;

        public:

            /**
             * @brief Construct a new J2 object.
             * 
             * @param[in] mu Central body gravitational parameter.
             * @param[in] J2 Central body J2-term.
             * @param[in] radius Central body radius.
             */
            J2(T mu, T J2, T radius);

            /**
             * @brief Calculate perturbing acceleration resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::array<T, 3> Total perturbing acceleration due to the J2-term.
             */
            std::array<T, 3> acceleration_total(T t, std::array<T, 3> R, std::array<T, 3> V) const override;

            /**
             * @brief Calculate perturbing acceleration resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<T> Total perturbing acceleration due to the J2-term.
             */
            std::vector<T> acceleration_total(T t, std::vector<T> R, std::vector<T> V) const override;

            /**
             * @brief Calculate perturbing potential resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return T Perturbing potential due to the J2-term.
             */
            T potential(T t, std::array<T, 3> R) const override;

            /**
             * @brief Calculate perturbing potential resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return T Perturbing potential due to the J2-term.
             */
            T potential(T t, std::vector<T> R) const override;

    };

}

#endif