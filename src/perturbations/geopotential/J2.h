#ifndef THAMES_PERTURBATIONS_GEOPOTENTIAL_J2
#define THAMES_PERTURBATIONS_GEOPOTENTIAL_J2

#include <array>
#include <vector>

#include "../baseperturbation.h"

using namespace thames::perturbations::baseperturbation;

namespace thames::perturbations::geopotential{

    ///////////
    // Reals //
    ///////////

    /**
     * @brief Class for the perturbation resulting from the J2-term.
     * 
     * @tparam T Numeric type.
     */
    template<class T>
    class J2 : public BasePerturbation<T> {
        private:

            /// Central body gravitational parameter.
            const T m_mu;

            /// Central body J2-term.
            const T m_J2;

            /// Central body radius.
            const T m_radius;

        public:

            /**
             * @brief Construct a new J2 object.
             * 
             * @param[in] mu Central body gravitational parameter.
             * @param[in] J2 Central body J2-term.
             * @param[in] radius Central body radius.
             */
            J2(const T& mu, const T& J2, const T& radius);

            /**
             * @brief Destroy the J2 object.
             * 
             */
            ~J2();

            ////////////
            // Arrays //
            ////////////

            /**
             * @brief Calculate perturbing acceleration resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::array<T, 3> Total perturbing acceleration due to the J2-term.
             */
            std::array<T, 3> acceleration_total(const T& t, const std::array<T, 3>& R, const std::array<T, 3>& V) const override;

            /**
             * @brief Calculate perturbing potential resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return T Perturbing potential due to the J2-term.
             */
            T potential(const T& t, const std::array<T, 3>& R) const override;

            /////////////
            // Vectors //
            /////////////

            /**
             * @brief Calculate perturbing acceleration resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<T> Total perturbing acceleration due to the J2-term.
             */
            std::vector<T> acceleration_total(const T& t, const std::vector<T>& R, const std::vector<T>& V) const override;

            /**
             * @brief Calculate perturbing potential resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return T Perturbing potential due to the J2-term.
             */
            T potential(const T& t, const std::vector<T>& R) const override;

    };

    /////////////////
    // Polynomials //
    /////////////////

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Class for the perturbation resulting from the J2-term.
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     */
    template<class T, template<class> class P>
    class J2Polynomial : public BasePerturbationPolynomial<T, P> {
        private:

            /// Central body gravitational parameter
            const T m_mu;

            /// Central body J2-term
            const T m_J2;

            /// Central body radius
            const T m_radius;

        public:

            /**
             * @brief Construct a new J2Polynomial object.
             * 
             * @param[in] mu Central body gravitational parameter.
             * @param[in] J2 Central body J2-term.
             * @param[in] radius Central body radius.
             */
            J2Polynomial(const T& mu, const T& J2, const T& radius);

            /**
             * @brief Destroy the J2Polynomial object
             * 
             */
            ~J2Polynomial();

            /**
             * @brief Calculate perturbing acceleration resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return std::vector<P<T>> Total perturbing acceleration due to the J2-term.
             */
            std::vector<P<T>> acceleration_total(const T& t, const std::vector<P<T>>& R, const std::vector<P<T>>& V) const override;

            /**
             * @brief Calculate perturbing potential resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return P<T> Perturbing potential due to the J2-term.
             */
            P<T> potential(const T& t, const std::vector<P<T>>& R) const override;

    };

    #endif

}

#endif