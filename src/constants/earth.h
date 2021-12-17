#ifndef THAMES_CONSTANTS_EARTH
#define THAMES_CONSTANTS_EARTH

namespace thames::constants::earth{

    /// Radius of the Earth [km]
    const double radius = 6378.0;

    /// Earth's gravitational parameter [km^3/s^2]
    const double mu = 3.986004414498200E+05;

    /// Earth's J2 term [-]
    const double J2 = 1.082635854E-03;

}

#endif