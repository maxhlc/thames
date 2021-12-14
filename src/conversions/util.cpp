#include "util.h"
#include "../types.h"

using namespace thames::types;

namespace thames::conversions::util {

    Matrix33 rot_x(const double &ang){
        // TODO: documentation

        // Declare rotation matrix
        Matrix33 rot;

        // Populate rotation matrix
        rot << 1.0,       0.0,      0.0,
               0.0,  cos(ang), sin(ang),
               0.0, -sin(ang), cos(ang);

        // Return rotation matrix
        return rot;
    }

    Matrix33 rot_y(const double &ang){
        // TODO: documentation

        // Declare rotation matrix
        Matrix33 rot;

        // Populate rotation matrix
        rot << cos(ang), 0.0, -sin(ang),
                    0.0, 1.0,       0.0,
               sin(ang), 0.0,  cos(ang);

        // Return rotation matrix
        return rot;
    }

    Matrix33 rot_z(const double &ang){
        // TODO: documentation

        // Declare rotation matrix
        Matrix33 rot;

        // Populate rotation matrix
        rot <<  cos(ang), sin(ang), 0.0,
               -sin(ang), cos(ang), 0.0,
                     0.0,      0.0, 1.0;

        // Return rotation matrix
        return rot;
    }

}