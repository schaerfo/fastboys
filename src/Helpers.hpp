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

#ifndef HELPERS_HPP
#define HELPERS_HPP

#include "Molecule.hpp"

#include <numbers>

#include <Eigen/Dense>

namespace helpers {

inline bool is_s_orbital(BasisFunction::Type t) {
    return t == BasisFunction::Type::s;
}

inline bool is_p_orbital(BasisFunction::Type t) {
    return t == BasisFunction::Type::px || t == BasisFunction::Type::py || t == BasisFunction::Type::pz;
}

inline double gaussian_product(double alpha, const Eigen::Vector3d& a, double beta, const Eigen::Vector3d& b) {
    return std::exp(- alpha * beta / (alpha + beta) * (a - b).squaredNorm());
}

inline double vector_component(const Eigen::Vector3d& vec, BasisFunction::Type comp) {
    assert(is_p_orbital(comp));
    auto i = static_cast<int>(comp) - 1;
    return vec[i];
}

inline double component_difference(const Eigen::Vector3d& a, const Eigen::Vector3d& b, BasisFunction::Type comp) {
    return vector_component(a, comp) - vector_component(b, comp);
}

inline double s1(double t) {
    return std::abs(t) < 1e-15 ? 2 / std::sqrt(std::numbers::pi) : std::erf(t) / t;
}

inline double s2(double t) {
    return std::abs(t) < 1e-5 ? -4.0 / (3 * std::sqrt(std::numbers::pi)) : (2 * std::pow(std::numbers::pi, -0.5) * t * std::exp(-t * t) - std::erf(t)) / (t * t * t);
}

inline double s3(double t) {
    return std::abs(t) < 1e-3 ? 8.0 / (5 * std::sqrt(std::numbers::pi)) : (3 * std::erf(t) - 2 * (3 * t + 2 * t * t * t) * std::exp(-t * t) / std::sqrt(std::numbers::pi)) / (t * t * t * t * t);
}

}

#endif //HELPERS_HPP
