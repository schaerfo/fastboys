#include <fstream>
#include <iostream>

#include "Molecule.hpp"

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
    return 0;
}
