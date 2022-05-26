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

#include <fstream>
#include <iomanip>
#include <string>

#include <nlohmann/json.hpp>

#include "../../include/io/json.h"

namespace thames::io::json {

    template<class T>
    void load(const std::string& filepath, Parameters<T>& parameters) {
        /// @todo Input checking

        // Open file stream
        std::ifstream filestream(filepath);

        // Load JSON
        nlohmann::json j;
        filestream >> j;

        // Load parameters
        parameters = j.get<Parameters<T>>();
    }
    template void load(const std::string&, Parameters<double>&);

    template<class T>
    void save(const std::string& filepath, const Parameters<T>& parameters) {
        // Open file stream
        std::ofstream filestream(filepath);

        // Construct JSON object
        nlohmann::json j = parameters;

        // Output JSON object
        filestream << std::setw(4) << j;
    }
    template void save(const std::string&, const Parameters<double>&); 

}