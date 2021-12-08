#include <array>
#include <iostream>

#include <eigen3/Eigen/Core>

#include "conversions/state.h"
#include "types.h"

using namespace thames::types;

int main(){
    Vector3 R = {7100.0f, 0.0f, 1300.0f};
    Vector3 V = {0.0f, 7.35f, 1.0f};
    double mu = 3.986e5;

    Vector6 keplerian = thames::conversions::state::cartesian_to_keplerian(R, V, mu);

    std::cout << keplerian[0] << " "
              << keplerian[1] << " "
              << keplerian[2] << " "
              << keplerian[3] << " "
              << keplerian[4] << " "
              << keplerian[5] << std::endl;

    std::cout << "THAMES executed successfully" << std::endl;

    return 0;
}