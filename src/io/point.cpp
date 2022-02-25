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
#include <iostream>
#include <vector>
#include <string>

#include "point.h"

namespace thames::io::point {

    template<class T>
    void load(const std::string filepath, T& tstart, T& tend, int& scid, int& statetype, std::vector<std::vector<T>>& states){
        // Clear states vector
        states.clear();

        // Open filestream
        std::ifstream filestream(filepath);

        // Throw error if unable to open file
        if(!filestream.is_open())
            throw std::invalid_argument("Could not open file: " + filepath);

        // Declare temporary variables for processing
        std::string line, word;
        std::vector<T> state;

        // Load header
        std::getline(filestream, line);
        std::stringstream str(line);

        // Import start time
        std::getline(str, word, ',');
        tstart = std::stod(word);

        // Import end time
        std::getline(str, word, ',');
        tend = std::stod(word);

        // Import spacecraft ID
        std::getline(str, word, ',');
        scid = std::stoi(word);

        // Import statetype
        std::getline(str, word, ',');
        statetype = std::stoi(word);

        // Import states
        while(std::getline(filestream, line)){
            // Clear temporary state vector
            state.clear();

            // Create string stream for current line
            std::stringstream str(line);

            // Extract each state variable and store in temporary state vector
            while(std::getline(str, word, ','))
                state.push_back(std::stod(word));

            // Add temporary state vector to overal states vector
            states.push_back(state);        
        }
    }
    template void load(const std::string, double&, double&, int&, int&, std::vector<std::vector<double>>&);

    template<class T>
    void load(const std::string filepath, T& tstart, T& tend, int& scid, std::vector<std::vector<T>>& states){
        int statetype;
        load(filepath, tstart, tend, scid, statetype, states);
    }
    template void load(const std::string, double&, double&, int&, std::vector<std::vector<double>>&);

    template<class T>
    void save(const std::string filepath, const T& tstart, const T& tend, const int& scid, const std::vector<std::vector<T>>& states, const unsigned int precision){
        // Open filestream
        std::ofstream filestream(filepath);

        // Throw error if unable to open file
        if(!filestream.is_open())
            throw std::invalid_argument("Could not open file: " + filepath);

        // Set precision
        filestream << std::setprecision(precision);

        // Print states
        for(std::size_t ii=0; ii<states.size(); ii++)
            filestream << states[ii][0] << ", " << states[ii][1] << ", " << states[ii][2] << ", " << states[ii][3] << ", " << states[ii][4] << ", " << states[ii][5] << "\n";

        // Print footer
        filestream << tstart << ", " << tend << ", " << scid;
    }
    template void save(const std::string, const double&, const double&, const int&, const std::vector<std::vector<double>>&, const unsigned int);

}