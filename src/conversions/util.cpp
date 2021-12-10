#include "util.h"
#include "../types.h"

using namespace thames::types;

namespace thames::conversions::util {

    Matrix33 rot_x(const double &ang){
        // TODO: documentation

        // Declare rotation matrix
        Matrix33 rot;

        // Populate rotation matrix
        rot(0, 0) = 1.0f;   rot(0, 1) = 0.0f;       rot(0, 2) = 0.0f;
        rot(1, 0) = 0.0f;   rot(1, 1) = cos(ang);   rot(1, 2) = sin(ang);
        rot(2, 0) = 0.0f;   rot(2, 1) = -sin(ang);  rot(2, 2) = cos(ang);

        // Return rotation matrix
        return rot;
    }

    Matrix33 rot_y(const double &ang){
        // TODO: documentation

        // Declare rotation matrix
        Matrix33 rot;

        // Populate rotation matrix
        rot(0, 0) = cos(ang);   rot(0, 1) = 0.0f;   rot(0, 2) = -sin(ang);
        rot(1, 0) = 0.0f;       rot(1, 1) = 1.0f;   rot(1, 2) = 0.0f;
        rot(2, 0) = sin(ang);   rot(2, 1) = 0.0f;   rot(2, 2) = cos(ang);

        // Return rotation matrix
        return rot;
    }

    Matrix33 rot_z(const double &ang){
        // TODO: documentation

        // Declare rotation matrix
        Matrix33 rot;

        // Populate rotation matrix
        rot(0, 0) = cos(ang);   rot(0, 1) = sin(ang);   rot(0, 2) = 0.0f;
        rot(1, 0) = -sin(ang);  rot(1, 1) = cos(ang);   rot(1, 2) = 0.0f;
        rot(2, 0) = 0.0f;       rot(2, 1) = 0.0f;       rot(2, 2) = 1.0f;

        // Return rotation matrix
        return rot;
    }

}