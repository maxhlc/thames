#include <vector>

#include "sampling.h"

namespace thames::util::sampling {

    template<class T>
    std::vector<T> linspace(const T a, const T b, const unsigned int n) {
        // Declare output vector
        std::vector<T> points;

        if(n == 1){
            // Return midpoint if only one point requested
            points.push_back((a + b)/2);
        } else {
            // Calculate step size
            T step = (b - a)/(n - 1);

            // Generate points
            for(unsigned int ii=0; ii<n; ii++)
                points.push_back(a + step*ii);
        }

        // Return points    
        return points;
    }
    template std::vector<double> linspace(const double, const double, const unsigned int);

    template<class T>
    std::vector<std::vector<T>> cartesian_permutations(const std::vector<std::vector<T>>& points) {
        // Declare output vector
        std::vector<std::vector<T>> permuations;

        // Iterate through provided points
        // TODO: fix horrible implementation to avoid so many nested loops
        for(std::size_t ix=0; ix < points[0].size(); ix++){
            for(std::size_t iy=0; iy < points[1].size(); iy++){
                for(std::size_t iz=0; iz < points[2].size(); iz++){
                    for(std::size_t ivx=0; ivx < points[3].size(); ivx++){
                        for(std::size_t ivy=0; ivy < points[4].size(); ivy++){
                            for(std::size_t ivz=0; ivz < points[5].size(); ivz++){
                                // Generate current permutation
                                std::vector<T> permutation = {
                                    points[0][ix],
                                    points[1][iy],
                                    points[2][iz],
                                    points[3][ivx],
                                    points[4][ivy],
                                    points[5][ivz]
                                };

                                // Store permuation
                                permuations.push_back(permutation);
                            }
                        }
                    }
                }
            }
        }

        // Return permutations
        return permuations;
    }
    template std::vector<std::vector<double>> cartesian_permutations(const std::vector<std::vector<double>>&);

}