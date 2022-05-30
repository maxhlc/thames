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

#ifndef THAMES_IO_JSON
#define THAMES_IO_JSON

#include "../settings/settings.h"

namespace thames::io::json {

    /**
     * @brief Load parameters from JSON
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-30
     * 
     * @tparam T Numeric type
     * @param[in] filepath Input file path
     * @param[out] parameters Parameters
     */
    template<class T>
    void load(const std::string& filepath, thames::settings::Parameters<T>& parameters);

    /**
     * @brief Save parameters to JSON
     * 
     * @author Max Hallgarten La Casta
     * @date 2022-05-30
     * 
     * @tparam T Numeric type
     * @param[in] filepath Input file path
     * @param[in] parameters Parameters
     */
    template<class T>
    void save(const std::string& filepath, const thames::settings::Parameters<T>& parameters);

}

#endif