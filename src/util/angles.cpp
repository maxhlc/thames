#include <math.h>

namespace thames::util::angles{

    double angle_wrap(double theta){
        // TODO: documentation
        // https://stackoverflow.com/questions/11498169/dealing-with-angle-wrap-in-c-code

        theta = fmod(theta + M_PI, 2*M_PI);

        if (theta < 0.0){
            theta += 2*M_PI;
        }

        return theta - M_PI;
    }

}