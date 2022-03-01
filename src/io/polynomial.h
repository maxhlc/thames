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

#ifndef THAMES_IO_POLYNOMIAL
#define THAMES_IO_POLYNOMIAL

#include <vector>
#include <string>

#include "../constants/statetypes.h"

namespace thames::io::polynomial {

    #ifdef THAMES_USE_SMARTUQ

    /**
     * @brief Load polynomial coefficients.
     * 
     * @todo This function is templated, however the conversion from string to a floating point representation currently only converts to double.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-03-01
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] filepath Filepath to be read.
     * @param[out] tstart Start time.
     * @param[out] tend End time.
     * @param[out] scid Spacecraft identifier.
     * @param[out] statetype State type (e.g. Cartesian, GEqOE, etc.)
     * @param[out] degree Degree of the polynomials. 
     * @param[out] polynomials Polynomial states.
     */
    template<class T, template<class> class P>
    void load(const std::string filepath, T& tstart, T& tend, int& scid, thames::constants::statetypes::StateTypes& statetype, int& degree, std::vector<P<T>>& polynomials);

    /**
     * @brief Save polynomial coefficients.
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-03-01
     * 
     * @tparam T Numeric type.
     * @tparam P Polynomial type.
     * @param[in] filepath Filepath to be read.
     * @param[in] tstart Start time.
     * @param[in] tend End time.
     * @param[in] scid Spacecraft identifier.
     * @param[in] statetype State type (e.g. Cartesian, GEqOE, etc.)
     * @param[in] degree Degree of the polynomials. 
     * @param[in] polynomials Polynomial states.
     * @param[in] precision Output precision.
     */
    template<class T, template<class> class P>
    void save(const std::string filepath, const T& tstart, const T& tend, const int& scid, const thames::constants::statetypes::StateTypes& statetype, const int& degree, const std::vector<P<T>>& polynomials, const unsigned int precision = 15);

    #endif

}

#endif