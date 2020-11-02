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

#include "Molecule.hpp"
#include "OneElectronIntegrals.hpp"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: fastboys <input.xyz> <basisset.json>\n";
        return 1;
    }

    std::ifstream input(argv[1]);
    if (!input.is_open()) {
        std::cerr << "Could not open file: " << argv[1] << '\n';
        return 2;
    }

    Molecule m(input);

    std::ifstream basis(argv[2]);
    if (!basis.is_open()) {
        std::cerr << "Could not open file: " << argv[2] << '\n';
        return 2;
    }
    auto b = m.construct_basis_set(basis);
    auto s = overlap(b);
    auto k = kinetic_energy(b);
    auto v = potential_energy(b, m);
    std::cout << s << "\n\n" << k << "\n\n" << v << '\n';
    return 0;
}
