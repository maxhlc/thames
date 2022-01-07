#ifndef THAMES_UTIL_ANGLES
#define THAMES_UTIL_ANGLES

namespace thames::util::angles{

    /**
     * @brief Wrap angle between -pi and pi.
     * 
     * @details Method described by: https://stackoverflow.com/a/11498248
     * 
     * @param[in] theta Angle.
     * @return double Wrapped angle.
     */
    double angle_wrap(double theta);

}

#endif