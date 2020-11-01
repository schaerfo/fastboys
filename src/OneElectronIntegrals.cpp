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

#include "OneElectronIntegrals.hpp"

#include "Helpers.hpp"

#include <numbers>

double add_primitive_integrals(const BasisFunction& r, const BasisFunction& s, auto func, auto obj) {
    double res = 0;
    for (const auto& [alpha, a]: r.primitives) {
        for (const auto& [beta, b]: s.primitives) {
            res += a * b * std::invoke(func, obj, alpha, r.center, r.type, beta, s.center, s.type);
        }
    }
    return res;
}

template<typename FunctionSet>
auto get_primitive_function(BasisFunction::Type r_type, BasisFunction::Type s_type) {
    if (helpers::is_s_orbital(r_type))
        return &FunctionSet::prim_ss;
    else if (helpers::is_s_orbital(s_type))
        return &FunctionSet::prim_ps;
    else
        return &FunctionSet::prim_pp;
}

template <typename FunctionSet>
Eigen::MatrixXd one_electron_matrix(const BasisSet& set, FunctionSet functions) {
    auto n = set.size();
    Eigen::MatrixXd res(n, n);
    res.setZero();
    for(std::size_t i=0; i<n; ++i) {
        for(std::size_t j=i; j<n; ++j) {
            const auto* r = &set[i];
            const auto* s = &set[j];

            if (helpers::is_s_orbital(r->type) && helpers::is_p_orbital(s->type))
                std::swap(r, s);

            res(i, j) = res(j, i) = add_primitive_integrals(*r, *s, get_primitive_function<FunctionSet>(r->type, s->type), functions);
        }
    }
    return res;
}

struct OverlapIntegrals {
    double prim_ss(double alpha, const Eigen::Vector3d& a, BasisFunction::Type,
                   double beta, const Eigen::Vector3d& b, BasisFunction::Type) {
        return helpers::gaussian_product(alpha, a, beta, b) * std::pow(std::numbers::pi / (alpha + beta), 1.5);
    }

    double prim_ps(double alpha, const Eigen::Vector3d& a, BasisFunction::Type i,
                   double beta, const Eigen::Vector3d& b, BasisFunction::Type) {
        return -helpers::gaussian_product(alpha, a, beta, b) * beta * std::pow(std::numbers::pi, 1.5) / std::pow(alpha + beta, 2.5) * helpers::component_difference(a, b, i);
    }

    double prim_pp(double alpha, const Eigen::Vector3d& a, BasisFunction::Type i,
                   double beta, const Eigen::Vector3d& b, BasisFunction::Type j) {
        return helpers::gaussian_product(alpha, a, beta, b) * std::pow(std::numbers::pi, 1.5) / std::pow(alpha + beta, 2.5) * (0.5 * (i == j) - alpha * beta / (alpha + beta) * helpers::component_difference(a, b, i) * helpers::component_difference(a, b, j));
    }
};

Eigen::MatrixXd overlap(const BasisSet& basis_set) {
    return one_electron_matrix(basis_set, OverlapIntegrals());
}
