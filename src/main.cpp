/**
 * Copyright (c) 2020 Christian Schärf
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <fstream>
#include <iostream>

#include <boost/program_options.hpp>
#include <omp.h>

#include "Molecule.hpp"
#include "OneElectronIntegrals.hpp"

namespace po = boost::program_options;

int main(int argc, char** argv) {
    po::options_description desc("fastboys options");
    desc.add_options()
            ("input", po::value<std::string>()->required(), "xyz input file (required)")
            ("basisset", po::value<std::string>()->required(), "JSON basisset file (required)");

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);
    } catch (const po::error& e) {
        std::cerr << e.what() << "\n\n" << desc << '\n';
        return 1;
    }

    std::ifstream input(vm["input"].as<std::string>());
    if (!input.is_open()) {
        std::cerr << "Could not open file: " << vm["input"].as<std::string>() << '\n';
        return 2;
    }

    Molecule m(input);

    std::ifstream basis(vm["basisset"].as<std::string>());
    if (!basis.is_open()) {
        std::cerr << "Could not open file: " << vm["basisset"].as<std::string>() << '\n';
        return 2;
    }
    auto b = m.construct_basis_set(basis);
    auto s = overlap(b);
    auto k = kinetic_energy(b);
    auto v = potential_energy(b, m);
    std::cout << s << "\n\n" << k << "\n\n" << v << '\n';
    return 0;
}
