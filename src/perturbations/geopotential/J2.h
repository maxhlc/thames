#ifndef THAMES_PERTURBATIONS_GEOPOTENTIAL_J2
#define THAMES_PERTURBATIONS_GEOPOTENTIAL_J2

#include "../baseperturbation.h"
#include "../../types.h"

using namespace thames::types;
using namespace thames::perturbations::baseperturbation;

namespace thames::perturbations::geopotential{

    /**
     * @brief Class for the perturbation resulting from the J2-term.
     * 
     * @tparam real Type for real numbers (e.g. float, double, etc.)
     * @tparam vector Type for vector (e.g. std::vector<double, 3>, Eigen::Vector3d)
     */
    template<class real, class vector>
    class J2 : public BasePerturbation<real, vector> {
        private:

            /// Central body gravitational parameter.
            real m_mu;

            /// Central body J2-term.
            real m_J2;

            /// Central body radius.
            real m_radius;

        public:

            /**
             * @brief Construct a new J2 object.
             * 
             * @param[in] mu Central body gravitational parameter.
             * @param[in] J2 Central body J2-term.
             * @param[in] radius Central body radius.
             */
            J2(real mu, real J2, real radius);

            /**
             * @brief Calculate perturbing acceleration resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @param[in] V Velocity vector.
             * @return vector Total perturbing acceleration due to the J2-term.
             */
            vector acceleration_total(real t, vector R, vector V) override;

            /**
             * @brief Calculate perturbing potential resulting from the J2-term. 
             * 
             * @param[in] t Current physical time.
             * @param[in] R Position vector.
             * @return real Perturbing potential due to the J2-term.
             */
            real potential(real t, vector R) override;

    };

}

#endif