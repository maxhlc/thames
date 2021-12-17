#include <math.h>

namespace thames::util::angles{

    double angle_wrap(double theta){

        theta = fmod(theta + M_PI, 2*M_PI);

        if (theta < 0.0){
            theta += 2*M_PI;
        }

        return theta - M_PI;
    }

}