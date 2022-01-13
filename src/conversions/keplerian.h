#ifndef THAMES_CONVERSIONS_KEPLERIAN
#define THAMES_CONVERSIONS_KEPLERIAN

#include "../types.h"

using namespace thames::types;

namespace thames::conversions::keplerian{

    /**
     * @brief Convert from Cartesian state to traditional Keplerian elements.
     * 
     * @param[in] RV Cartesian state.
     * @param[in] mu Gravitational parameter.
     * @return Vector6 Keplerian elements state.
     */
    template<class real, class vector3, class vector6>
    vector6 cartesian_to_keplerian(const vector6 &RV, const real &mu);

    /**
     * @brief Convert from traditional Keplerian elements to Cartesian state.
     * 
     * @param[in] keplerian Keplerian elements state.
     * @param[in] mu Gravitational parameter.
     * @return Vector6 Cartesian state.
     */
    template<class real, class vector3, class vector6>
    vector6 keplerian_to_cartesian(const vector6 &keplerian, const real &mu);

}

#endif