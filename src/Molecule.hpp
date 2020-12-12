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
#include <nlohmann/json.hpp>

class Atom {
public:
    Atom(const std::string& symbol, Eigen::Vector3d pos);

    Eigen::Vector3d pos_;
    unsigned charge_;

    bool operator==(const Atom& other) const;
};

struct BasisFunction {
    struct PrimitiveFunction {
        double exponent, coefficient;
        bool operator==(const PrimitiveFunction& other) const = default;
    };
    enum class Type {
        s = 0b0, px = 0b1, py = 0b10, pz = 0b11
    };
    using PrimitiveVec = boost::container::small_vector<PrimitiveFunction, 6>;

    Eigen::Vector3d center;
    PrimitiveVec primitives;
    Type type;

    bool operator==(const BasisFunction& other) const = default;
};

std::ostream& operator<<(std::ostream& os, const BasisFunction& function);

using BasisSet = std::vector<BasisFunction>;
std::ostream& operator<<(std::ostream& os, const BasisSet& basis_set);

class Molecule {
public:
    Molecule() = default;
    explicit Molecule(std::vector<Atom> atoms);

    explicit Molecule(std::istream& xyz);

    BasisSet construct_basis_set(std::istream& basisset_json) const;
    BasisSet construct_basis_set(const nlohmann::json& basisset) const;

    const std::vector<Atom>& get_atoms() const {
        return atoms_;
    }

    unsigned int get_occupied_orbitals() const {
        return occupied_orbitals_;
    }

    double electronic_energy(const Eigen::MatrixXd& h_mo, const Eigen::MatrixXd& f_mo) const;
    double nuclear_energy() const;

    bool operator==(const Molecule& other) const = default;

private:
    std::vector<Atom> atoms_;
    unsigned occupied_orbitals_;

    unsigned calc_occupied_orbitals() const;
};

std::ostream& operator<<(std::ostream& os, const Atom& atom);
std::ostream& operator<<(std::ostream& os, const Molecule& molecule);

#endif //MOLECULE_HPP
