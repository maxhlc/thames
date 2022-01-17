#ifndef THAMES_UTIL_ANGLES
#define THAMES_UTIL_ANGLES

namespace thames::util::angles{

    /**
     * @brief Wrap angle between -pi and pi.
     * 
     * @details Method described by: https://stackoverflow.com/a/11498248
     * 
     * @tparam T Numeric type.
     * @param[in] theta Angle.
     * @return T Wrapped angle.
     */
    template<class T>
    T angle_wrap(T theta);

}

#endif