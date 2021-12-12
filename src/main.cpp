#include <array>
#include <iostream>

#include <Eigen/Core>

#include "conversions/state.h"
#include "types.h"

using namespace thames::types;

double U(double t, Vector3 R){
    return 0.0f;
}

int main(){
    Vector6 RV;
    RV << 7100.0, 0.0, 1300.0, 0.3, 7.35, 1.0;
    double mu = 3.986e5;

    Vector3 R = RV(Eigen::seq(0,2));
    Vector3 V = RV(Eigen::seq(3,5));

    Vector6 keplerian = thames::conversions::state::cartesian_to_keplerian(RV, mu);
    Vector6 state = thames::conversions::state::keplerian_to_cartesian(keplerian, mu);

    Vector3 R2 = state(Eigen::seq(0,2));
    Vector3 V2 = state(Eigen::seq(3,5));

    double pos_error = (R-R2).norm();
    double vel_error = (V-V2).norm();

    double t = 0.0;

    Vector6 geqoe = thames::conversions::state::cartesian_to_geqoe(t, RV, mu, U);
    Vector6 state2 = thames::conversions::state::geqoe_to_cartesian(t, geqoe, mu, U);

    std::cout << geqoe << std::endl;

    Vector3 R3 = state2(Eigen::seq(0,2));
    Vector3 V3 = state2(Eigen::seq(3,5));

    double pos2_error = (R-R3).norm();
    double vel2_error = (V-V3).norm();

    std::cout << "Keplerian pos/vel error: " << pos_error << " " << vel_error << std::endl;
    std::cout << "GEqOE pos/vel error: " << pos2_error << " " << vel2_error << std::endl;
    std::cout << "THAMES executed successfully" << std::endl;

    return 0;
}