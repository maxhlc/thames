#include <array>
#include <iostream>

#include <Eigen/Core>

#include "conversions/state.h"
#include "types.h"

using namespace thames::types;

int main(){
    Vector6 RV;
    RV << 7100.0f, 0.0f, 1300.0f, 0.0f, 7.35f, 1.0f;
    double mu = 3.986e5;

    Vector6 keplerian = thames::conversions::state::cartesian_to_keplerian(RV, mu);
    Vector6 state = thames::conversions::state::keplerian_to_cartesian(keplerian, mu);

    std::cout << keplerian << std::endl;
    std::cout << state << std::endl;

    std::cout << "THAMES executed successfully" << std::endl;

    return 0;
}