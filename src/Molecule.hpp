#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include <istream>
#include <vector>
#include <Eigen/Dense>

class Atom {
public:
    Atom(std::string_view symbol, Eigen::Vector3d pos);

    std::string symbol_;
    Eigen::Vector3d pos_;
    unsigned charge_;
};

class Molecule {
public:
    explicit Molecule(std::istream& xyz);

private:
    std::vector<Atom> atoms_;
};

#endif //MOLECULE_HPP
