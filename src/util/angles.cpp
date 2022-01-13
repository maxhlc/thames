#include <cmath>

#include "angles.h"

namespace thames::util::angles{

    template<class real>
    real angle_wrap(real theta){

        theta = fmod(theta + M_PI, 2*M_PI);

        if (theta < 0.0){
            theta += 2*M_PI;
        }

        return theta - M_PI;
    }
    template double angle_wrap<double>(double);

}