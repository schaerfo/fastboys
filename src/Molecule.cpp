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

#include "Molecule.hpp"

#include <bitset>
#include <numbers>
#include <unordered_map>

#include <fmt/core.h>

constexpr double ang2bohr = 1.89;

const std::unordered_map<std::string, unsigned> atom_charges = {{"H" , 1},
                                                                {"He", 2},
                                                                {"Li", 3},
                                                                {"Be", 4},
                                                                {"B" , 5},
                                                                {"C" , 6},
                                                                {"N" , 7},
                                                                {"O" , 8},
                                                                {"F" , 9},
                                                                {"Ne", 10},
                                                                {"Na", 11},
                                                                {"Mg", 12},
                                                                {"Al", 13},
                                                                {"Si", 14},
                                                                {"P" , 15},
                                                                {"S" , 16},
                                                                {"Cl", 17},
                                                                {"Ar", 18},
                                                                {"Sc", 21}};

Atom::Atom(const std::string& symbol, Eigen::Vector3d pos):
    pos_(std::move(pos))
{
    try {
        charge_ = atom_charges.at(symbol);
    } catch (std::out_of_range &) {
        throw std::runtime_error(fmt::format("Unknown element: {}", symbol));
    }
}

bool Atom::operator==(const Atom& other) const {
    return other.charge_ == charge_ && other.pos_.isApprox(pos_, 1e-10);
}

bool BasisFunction::PrimitiveFunction::operator==(const PrimitiveFunction& other) const {
    return std::abs(exponent - other.exponent) < 1e-12 && std::abs(coefficient - other.coefficient) < 1e-12;
}

Molecule::Molecule(std::vector<Atom> atoms)
        : atoms_(std::move(std::move(atoms)))
        , occupied_orbitals_(calc_occupied_orbitals())
{}

Molecule::Molecule(std::istream& xyz) {
    int atom_count;
    xyz >> atom_count;
    atoms_.reserve(atom_count);
    xyz.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    xyz.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    for (int i=0; i<atom_count; ++i) {
        std::string symbol;
        double x, y, z;
        xyz >> symbol >> x >> y >> z;
        atoms_.emplace_back(symbol, ang2bohr * Eigen::Vector3d(x, y, z));
    }

    occupied_orbitals_ = calc_occupied_orbitals();
}

double norm_s(double exponent) {
    return std::pow(2 * exponent / std::numbers::pi, 0.75);
}

double norm_p(double exponent) {
    return std::pow(128 * std::pow(exponent, 5) / std::pow(std::numbers::pi, 3), 0.25);
}

void norm_primitives(BasisFunction::PrimitiveVec& primitives, int l) {
    assert(l == 0 || l == 1);
    auto norm_func = l == 0 ? norm_s : norm_p;
    for (auto& [curr_exp, curr_coeff]: primitives)
        curr_coeff *= norm_func(curr_exp);
}

template <typename OutputIt>
void add_basisfunctions(OutputIt out, const nlohmann::json& shells, const Eigen::Vector3d& center) {
    for (const auto& [foo, curr_shell]: shells.items()) {
        for (unsigned i=0; i<curr_shell["coefficients"].size(); ++i) {
            int curr_l = curr_shell["angular_momentum"][i];
            if (curr_l != 0 && curr_l != 1)
                throw std::runtime_error(fmt::format("Unsupported angular quantum number: {}", curr_l));

            auto curr_coefficients = curr_shell["coefficients"][i];
            BasisFunction::PrimitiveVec primitives;
            for (unsigned i=0; i<curr_coefficients.size(); ++i)
                primitives.emplace_back(BasisFunction::PrimitiveFunction{
                    std::stod(static_cast<const std::string&>(curr_shell["exponents"][i])),
                    std::stod(static_cast<const std::string&>(curr_coefficients[i]))});

            norm_primitives(primitives, curr_l);
            if (curr_l == 0)
                *out++ = BasisFunction{center, primitives, BasisFunction::Type::s};
            else {
                *out++ = BasisFunction{center, primitives, BasisFunction::Type::px};
                *out++ = BasisFunction{center, primitives, BasisFunction::Type::py};
                *out++ = BasisFunction{center, primitives, BasisFunction::Type::pz};
            }
        }
    }
}

BasisSet Molecule::construct_basis_set(std::istream& basisset_json) const {
    nlohmann::json basisset;
    basisset_json >> basisset;
    return construct_basis_set(basisset);
}

BasisSet Molecule::construct_basis_set(const nlohmann::json& basisset) const {
    BasisSet res;
    for (const auto& curr_atom: atoms_)
        add_basisfunctions(std::back_inserter(res), basisset["elements"][std::to_string(curr_atom.charge_)]["electron_shells"], curr_atom.pos_);
    return res;
}

double Molecule::electronic_energy(const Eigen::MatrixXd& h_mo, const Eigen::MatrixXd& f_mo) const{
    double result = 0;
    assert(h_mo.cols() == h_mo.rows() && h_mo.cols() == f_mo.rows() && h_mo.cols() == f_mo.cols());
    for (Eigen::Index i=0; i < occupied_orbitals_; ++i) {
        result += h_mo(i, i) + f_mo(i, i);
    }
    return result;
}

unsigned Molecule::calc_occupied_orbitals() const {
    return std::accumulate(atoms_.begin(), atoms_.end(), 0,
                           [](unsigned sum, const Atom& a){return sum + a.charge_;}) / 2;
}

double Molecule::nuclear_energy() const {
    double result = 0;
    for (std::size_t i=0; i<atoms_.size(); ++i) {
        for (std::size_t j=i+1; j<atoms_.size(); ++j) {
            result += atoms_[i].charge_ * atoms_[j].charge_ / (atoms_[i].pos_ - atoms_[j].pos_).norm();
        }
    }
    return result;
}

std::ostream& operator<<(std::ostream& os, const Atom& atom) {
    return os << "Charge: " << atom.charge_ << " Coordinates: " << atom.pos_[0] << ' ' << atom.pos_[1] << ' ' << atom.pos_[2];
}

std::ostream& operator<<(std::ostream& os, const Molecule& m) {
    os << "Number of electrons: " << m.get_occupied_orbitals() << '\n';
    for (auto atoms = m.get_atoms(); const auto& curr_atom: atoms) {
        os << curr_atom << '\n';
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const BasisFunction& function) {
    os << "Type: " << std::bitset<2>(static_cast<int>(function.type)) << "; Center: " << function.center.transpose() << "; Primitives: ";
    for (auto [exp, coeff]: function.primitives) {
        os << exp << ' ' << coeff << "; ";
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, const BasisSet& basis_set) {
    for (const auto& curr_func: basis_set) {
        os << curr_func << '\n';
    }
    return os;
}
