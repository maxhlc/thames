#ifndef THAMES_UTIL_SAMPLING
#define THAMES_UTIL_SAMPLING

namespace thames::util::sampling {

    /**
     * @brief Calculate evenly spaced points over the [a,b] interval.
     * 
     * Returns the mid-point if one point is requested.
     * 
     * @tparam T Numeric type.
     * @param[in] a Start value.
     * @param[in] b End value.
     * @param[in] n Number of points.
     * @return std::vector<T> Evenly space points.
     */
    template<class T>
    std::vector<T> linspace(const T a, const T b, const unsigned int n);

    /**
     * @brief Generate permuations of Cartesian states from sets of points in each state variable.
     * 
     * @tparam T Numeric type.
     * @param[in] points Vector of vectors of ranges for each Cartesian state (i.e. X, Y, Z, VX, VY, VZ)
     * @return std::vector<std::vector<T>> 
     */
    template<class T>
    std::vector<std::vector<T>> cartesian_permutations(const std::vector<std::vector<T>>& points);

}

#endif