/*
MIT License

Copyright (c) 2021-2022 Max Hallgarten La Casta

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

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