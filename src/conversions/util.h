#ifndef THAMES_CONVERSIONS_UTIL
#define THAMES_CONVERSIONS_UTIL

#include "../types.h"

using namespace thames::types;

namespace thames::conversions::util{

    /**
     * @brief Generate coordinate system rotation matrix about the x-axis.
     * 
     * @param[in] ang Rotation angle.
     * @return Matrix33 Rotation matrix.
     */
    Matrix33 rot_x(const double &ang);

    /**
     * @brief Generate coordinate system rotation matrix about the y-axis.
     * 
     * @param[in] ang Rotation angle.
     * @return Matrix33 Rotation matrix.
     */
    Matrix33 rot_y(const double &ang);

    /**
     * @brief Generate coordinate system rotation matrix about the z-axis.
     * 
     * @param[in] ang Rotation angle.
     * @return Matrix33 Rotation matrix.
     */
    Matrix33 rot_z(const double &ang);

}

#endif