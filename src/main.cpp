#include <array>
#include <iostream>

#include <Eigen/Core>

#include "types.h"
#include "conversions/dimensional.h"
#include "conversions/state.h"
#include "propagators/cowell.h"


using namespace thames::types;

Vector3 F_func(double t, Vector3 R, Vector3 V){
    Vector3 F;
    F << 0.0, 0.0, 0.0;
    return F;
}

double U(double t, Vector3 R){
    return 0.0f;
}

int main(){
    Vector6 RV;
    RV << 7100.0, 0.0, 1300.0, 0.0, 7.35, 1.0;
    double mu = 3.986004414498200E+05;
    double t = 0.0;

    Vector3 R = RV(Eigen::seq(0,2));
    Vector3 V = RV(Eigen::seq(3,5));

    Vector6 keplerian = thames::conversions::state::cartesian_to_keplerian(RV, mu);
    Vector6 state = thames::conversions::state::keplerian_to_cartesian(keplerian, mu);

    Vector6 geqoe = thames::conversions::state::cartesian_to_geqoe(t, RV, mu, U);
    Vector6 state2 = thames::conversions::state::geqoe_to_cartesian(t, geqoe, mu, U);

    Vector3 R2 = state(Eigen::seq(0,2));
    Vector3 V2 = state(Eigen::seq(3,5));

    double pos_error = (R-R2).norm();
    double vel_error = (V-V2).norm();

    Vector3 R3 = state2(Eigen::seq(0,2));
    Vector3 V3 = state2(Eigen::seq(3,5));

    double pos2_error = (R-R3).norm();
    double vel2_error = (V-V3).norm();

    std::cout << "Keplerian pos/vel error: " << pos_error << " " << vel_error << std::endl;
    std::cout << "GEqOE pos/vel error: " << pos2_error << " " << vel2_error << std::endl;

    Vector6 state_prop = thames::propagators::cowell::propagate(0.0, 86400.0, 1.0, state, mu, F_func);

    std::cout << state_prop << std::endl;

    std::cout << "THAMES executed successfully" << std::endl;

    return 0;
}