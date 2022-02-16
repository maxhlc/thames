#ifndef THAMES_IO_POINT
#define THAMES_IO_POINT

#include <vector>
#include <string>

namespace thames::io::point {

    /**
     * @brief Function to load Cartesian state vectors from a text file.
     * 
     * @todo This function is templated, however the conversion from string to a floating point representation currently only converts to double.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-16
     * 
     * @tparam T Numeric type.
     * @param[in] filepath Filepath to be read.
     * @param[out] tstart Start time.
     * @param[out] tend End time.
     * @param[out] scid Spacecraft identifier.
     * @param[out] states State vectors.
     */
    template<class T>
    void load(const std::string filepath, T& tstart, T& tend, int& scid, std::vector<std::vector<T>>& states);

    /**
     * @brief Function to save Cartesian state vectors to a text file.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-02-16
     * 
     * @tparam T Numeric type.
     * @param[in] filepath Filepath to be read.
     * @param[in] tstart Start time.
     * @param[in] tend End time.
     * @param[in] scid Spacecraft identifier.
     * @param[in] states State vectors.
     * @param[in] precision Output precision.
     */
    template<class T>
    void save(const std::string filepath, const T& tstart, const T& tend, const int& scid, const std::vector<std::vector<T>>& states, const unsigned int precision = 15);

}

#endif