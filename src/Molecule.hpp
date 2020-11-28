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

#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include <istream>
#include <vector>

#include <boost/container/small_vector.hpp>
#include <Eigen/Dense>

class Atom {
public:
    Atom(std::string_view symbol, Eigen::Vector3d pos);

    std::string symbol_;
    Eigen::Vector3d pos_;
    unsigned charge_;
};

struct BasisFunction {
    struct PrimitiveFunction {
        double exponent, coefficient;
    };
    enum class Type {
        s = 0b0, px = 0b1, py = 0b10, pz = 0b11
    };
    using PrimitiveVec = boost::container::small_vector<PrimitiveFunction, 6>;

    Eigen::Vector3d center;
    PrimitiveVec primitives;
    Type type;
};

using BasisSet = std::vector<BasisFunction>;

class Molecule {
public:
    Molecule() = default;
    Molecule(std::vector<Atom> atoms): atoms_(atoms) {}

    explicit Molecule(std::istream& xyz);

    BasisSet construct_basis_set(std::istream& basisset_json) const;

    const std::vector<Atom>& get_atoms() const {
        return atoms_;
    }

    void setAtoms(const std::vector<Atom>& atoms) {
        atoms_ = atoms;
    }

private:
    std::vector<Atom> atoms_;
};

#endif //MOLECULE_HPP
