#include "Molecule.hpp"

#include <unordered_map>
#include <fmt/core.h>

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
                                                                {"Ar", 18}};

Atom::Atom(std::string_view symbol, Eigen::Vector3d pos):
    symbol_(symbol),
    pos_(std::move(pos))
{
    try {
        charge_ = atom_charges.at(symbol_);
    } catch (std::out_of_range &) {
        throw std::runtime_error(fmt::format("Unknown element: {}", symbol));
    }
}

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
        atoms_.emplace_back(symbol, Eigen::Vector3d(x, y, z));

    }
}
