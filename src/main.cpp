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

#include <chrono>
#include <fstream>
#include <iostream>
#include <type_traits>

#include <boost/program_options.hpp>
#include <fmt/chrono.h>
#include <fmt/format.h>
#include <omp.h>

#include "DensityMatrix.hpp"
#include "FockMatrix.hpp"
#include "Molecule.hpp"
#include "OneElectronIntegrals.hpp"
#include "TwoElectronIntegral.hpp"

namespace po = boost::program_options;

constexpr double epsilon_e = 1e-9;
constexpr double epsilon_p = 1e-5;

int main(int argc, char** argv) {
    using namespace std::chrono;
    using MyClock = std::conditional_t<high_resolution_clock::is_steady, high_resolution_clock, steady_clock>;
    using MyMillisec = duration<double, std::milli>;

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

    auto start = MyClock::now();
    auto s = overlap(b);
    auto k = kinetic_energy(b);
    auto v = potential_energy(b, m);
    auto end = MyClock::now();
    fmt::print("Calculating one-electron integrals finished after {}\n", duration_cast<MyMillisec>(end-start));

    start = MyClock::now();
    auto two_electron_integrals = calculate_two_electron_integrals(b);
    end = MyClock::now();
    fmt::print("Calculating two-electron integrals finished after {}\n", duration_cast<MyMillisec>(end-start));

    auto c = initial_coefficients(s);
    DensityMatrix p(b.size(), m.get_occupied_orbitals());
    Eigen::MatrixXd h = k + v;
    Eigen::MatrixXd f(b.size(), b.size());
    f.setZero();
    auto f_mo = transform_matrix(f, c);
    auto old_energy = m.electronic_energy(transform_matrix(h, c), f_mo);

    auto step_count = 0u;
    start = MyClock::now();
    while (true) {
        ++step_count;
        p.updateDensity(c);
        f = h + electron_repulsion_matrix(two_electron_integrals, p);
        f_mo = transform_matrix(f, c);
        auto energy = m.electronic_energy(transform_matrix(h, c), f_mo);

        auto delta_e = std::abs(old_energy - energy);
        auto delta_p = p.difference_norm();

        fmt::print("ΔP = {}; ΔE = {}\n", delta_e, delta_p);
        if (delta_e < epsilon_e && delta_p < epsilon_p) {
            end = MyClock::now();
            fmt::print("HF converged after {} steps. Total energy: {}. Duration: {}", step_count, energy + m.nuclear_energy(), duration_cast<MyMillisec>(end-start));
            break;
        }

        c = update_coefficients(c, f_mo);
        old_energy = energy;
    }

    return 0;
}
