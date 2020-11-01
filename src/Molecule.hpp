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
    //using PrimitiveVec = boost::container::small_vector<PrimitiveFunction, 6>;
    using PrimitiveVec = std::vector<PrimitiveFunction>;

    Eigen::Vector3d center;
    PrimitiveVec primitives;
    Type type;
};

using BasisSet = std::vector<BasisFunction>;

class Molecule {
public:
    explicit Molecule(std::istream& xyz);

    BasisSet construct_basis_set(std::istream& basisset_json) const;

private:
    std::vector<Atom> atoms_;
};

#endif //MOLECULE_HPP
