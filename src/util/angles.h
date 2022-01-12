#ifndef THAMES_UTIL_ANGLES
#define THAMES_UTIL_ANGLES

namespace thames::util::angles{

    /**
     * @brief Wrap angle between -pi and pi.
     * 
     * @details Method described by: https://stackoverflow.com/a/11498248
     * 
     * @tparam real Type for a real number.
     * @param[in] theta Angle.
     * @return real Wrapped angle.
     */
    template<class real>
    real angle_wrap(real theta);

}

#endif