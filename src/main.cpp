#include <fstream>
#include <iostream>

#include "Molecule.hpp"

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: fastboys <input.xyz>\n";
        return 1;
    }

    std::ifstream f(argv[1]);
    if (!f.is_open()) {
        std::cerr << "Could not open file: " << argv[1] << '\n';
        return 2;
    }

    Molecule m(f);
    return 0;
}
