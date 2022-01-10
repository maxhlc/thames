#ifndef THAMES_PERTURBATIONS_GEOPOTENTIAL
#define THAMES_PERTURBATIONS_GEOPOTENTIAL

#include "baseperturbation.h"
#include "../types.h"

using namespace thames::types;
using namespace thames::perturbations::baseperturbation;

namespace thames::perturbations::geopotential{

    class J2 : public BasePerturbation {
        private:

            /// Central body gravitational parameter.
            double m_mu;

            /// Central body J2-term.
            double m_J2;

            /// Central body radius.
            double m_radius;

        public:

            /**
             * @brief Construct a new J2 object.
             * 
             * @param[in] mu Central body gravitational parameter.
             * @param[in] J2 Central body J2-term.
             * @param[in] radius Central body radius.
             */
            J2(double mu, double J2, double radius);

            /**
             * @brief Calculate perturbing acceleration resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return Vector3 Total perturbing acceleration due to the J2-term.
             */
            Vector3 acceleration_total(double t, Vector3 R, Vector3 V) override;

            /**
             * @brief Calculate perturbing potential resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return double Perturbing potential due to the J2-term.
             */
            double potential(double t, Vector3 R) override;

    };

}

#endif