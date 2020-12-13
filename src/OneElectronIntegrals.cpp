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
#include "OneElectronIntegralDetail.hpp"

#include "Helpers.hpp"

#include <numbers>

template <typename IntegralFunc>
double add_primitive_integrals(const BasisFunction& r, const BasisFunction& s, IntegralFunc func) {
    double res = 0;
    for (const auto& [alpha, a]: r.primitives) {
        for (const auto& [beta, b]: s.primitives) {
            res += a * b * func(alpha, beta);
        }
    }
    return res;
}

template<typename FunctionSet>
auto get_primitive_function(BasisFunction::Type r_type, BasisFunction::Type s_type) {
    if (helpers::is_s_orbital(r_type))
        return &FunctionSet::integral_ss;
    else if (helpers::is_s_orbital(s_type))
        return &FunctionSet::integral_ps;
    else
        return &FunctionSet::integral_pp;
}

void order_orbital_tuple(const BasisFunction*& r, const BasisFunction*& s) {
    if (helpers::is_s_orbital(r->type) && helpers::is_p_orbital(s->type))
        std::swap(r, s);
}

template <typename FunctionSet>
Eigen::MatrixXd one_electron_matrix(const BasisSet& set, FunctionSet functions) {
    auto n = set.size();
    Eigen::MatrixXd res(n, n);
    res.setZero();
#ifdef ENABLE_OPENMP
    #pragma omp parallel for schedule(dynamic, 1) default(none) shared(set, n, res, functions)
#endif
    for(Eigen::Index i = n-1; i >= 0; --i) {
        for(std::size_t j=i; j<n; ++j) {
            const auto* r = &set[i];
            const auto* s = &set[j];
            order_orbital_tuple(r, s);

            res(i, j) = res(j, i) = std::invoke(get_primitive_function<FunctionSet>(r->type, s->type), functions, *r, *s);
        }
    }
    return res;
}

double OverlapIntegrals::integral_ss(const BasisFunction& r, const BasisFunction& s) {
    assert(helpers::is_s_orbital(r.type) && helpers::is_s_orbital(s.type));
    return add_primitive_integrals(r, s, [a = r.center, b = s.center](auto alpha, auto beta) {
        return helpers::gaussian_product(alpha, a, beta, b) * std::pow(std::numbers::pi / (alpha + beta), 1.5);
    });
}

double OverlapIntegrals::integral_ps(const BasisFunction& r, const BasisFunction& s) {
    assert(helpers::is_p_orbital(r.type) && helpers::is_s_orbital(s.type));
    return add_primitive_integrals(r, s, [a = r.center, i = r.type, b = s.center](auto alpha, auto beta) {
        return -helpers::gaussian_product(alpha, a, beta, b) * beta * std::pow(std::numbers::pi, 1.5) /
               std::pow(alpha + beta, 2.5) * helpers::component_difference(a, b, i);
    });
}

double OverlapIntegrals::integral_pp(const BasisFunction& r, const BasisFunction& s) {
    assert(helpers::is_p_orbital(r.type) && helpers::is_p_orbital(s.type));
    return add_primitive_integrals(r, s, [a = r.center, i = r.type, b = s.center, j = s.type](auto alpha, auto beta) {
        return helpers::gaussian_product(alpha, a, beta, b) * std::pow(std::numbers::pi, 1.5) /
               std::pow(alpha + beta, 2.5) * (0.5 * (i == j) - alpha * beta / (alpha + beta) *
                                                               helpers::component_difference(a, b, i) *
                                                               helpers::component_difference(a, b, j));
    });
}

double KineticEnergyIntegrals::integral_ss(const BasisFunction& r, const BasisFunction& s) {
    assert(helpers::is_s_orbital(r.type) && helpers::is_s_orbital(s.type));
    return add_primitive_integrals(r, s, [a = r.center, b = s.center](auto alpha, auto beta) {
        return helpers::gaussian_product(alpha, a, beta, b) * alpha * beta * std::pow(std::numbers::pi, 1.5) /
               std::pow(alpha + beta, 2.5)
               * (3 - 2 * alpha * beta / (alpha + beta) * (a - b).squaredNorm());
    });
}

double KineticEnergyIntegrals::integral_ps(const BasisFunction& r, const BasisFunction& s) {
    assert(helpers::is_p_orbital(r.type) && helpers::is_s_orbital(s.type));
    return add_primitive_integrals(r, s, [a = r.center, i = r.type, b = s.center](auto alpha, auto beta) {
        return helpers::gaussian_product(alpha, a, beta, b) * alpha * beta * beta * std::pow(std::numbers::pi, 1.5) /
               std::pow(alpha + beta, 3.5)
               * (2 * alpha * beta / (alpha + beta) * (a - b).squaredNorm() - 5) *
               helpers::component_difference(a, b, i);
    });
}

double KineticEnergyIntegrals::integral_pp(const BasisFunction& r, const BasisFunction& s) {
    assert(helpers::is_p_orbital(r.type) && helpers::is_p_orbital(s.type));
    return add_primitive_integrals(r, s, [a = r.center, i = r.type, b = s.center, j = s.type](auto alpha, auto beta) {
        return helpers::gaussian_product(alpha, a, beta, b) * alpha * beta * std::pow(std::numbers::pi, 1.5) /
               std::pow(alpha + beta, 3.5)
               * ((i == j) * (2.5 - alpha * beta / (alpha + beta) * (a - b).squaredNorm())
                  +
                  alpha * beta / (alpha + beta) * (2 * alpha * beta / (alpha + beta) * (a - b).squaredNorm() - 7) *
                  helpers::component_difference(a, b, i) * helpers::component_difference(a, b, j));
    });
}

double NuclearPotentialIntegrals::integral_ss(const BasisFunction& r, const BasisFunction& s) {
    assert(helpers::is_s_orbital(r.type) && helpers::is_s_orbital(s.type));
    return add_primitive_integrals(r, s, [this, a = r.center, b = s.center](auto alpha, auto beta) {
        return nuclear_loop([&, g = helpers::gaussian_product(alpha, a, beta, b),
                                    eta = alpha + beta,
                                    r_raw = (alpha * a + beta * b) / (alpha + beta)](const Eigen::Vector3d& pos) {
            Eigen::Vector3d r = r_raw - pos;
            double t = std::sqrt(eta) * r.norm();
            return g * std::pow(std::numbers::pi, 1.5) / eta * helpers::s1(t);
        });
    });
}

double NuclearPotentialIntegrals::integral_ps(const BasisFunction& r, const BasisFunction& s) {
    assert(helpers::is_p_orbital(r.type) && helpers::is_s_orbital(s.type));
    return add_primitive_integrals(r, s, [this, a = r.center, i = r.type, b = s.center](auto alpha, auto beta) {
        return nuclear_loop([&, g = helpers::gaussian_product(alpha, a, beta, b),
                                    eta = alpha + beta,
                                    r_raw = (alpha * a + beta * b) / (alpha + beta)](const Eigen::Vector3d& pos) {
            Eigen::Vector3d r = r_raw - pos;
            double t = std::sqrt(eta) * r.norm();
            return g * std::pow(std::numbers::pi, 1.5) / (2 * eta) *
                   (helpers::vector_component(r, i) * helpers::s2(t) -
                    2 * beta * helpers::component_difference(a, b, i) / eta * helpers::s1(t));
        });
    });
}

double NuclearPotentialIntegrals::integral_pp(const BasisFunction& r, const BasisFunction& s) {
    assert(helpers::is_p_orbital(r.type) && helpers::is_p_orbital(s.type));
    return add_primitive_integrals(r, s, [this, a = r.center, i = r.type, b = s.center, j = s.type](auto alpha, auto beta) {
        return nuclear_loop([&, g = helpers::gaussian_product(alpha, a, beta, b),
                                    eta = alpha + beta,
                                    r_raw = (alpha * a + beta * b) / (alpha + beta)](const Eigen::Vector3d& pos) {
            Eigen::Vector3d r = r_raw - pos;
            double t = std::sqrt(eta) * r.norm();
            double r_i = helpers::vector_component(r, i), r_j = helpers::vector_component(r, j);
            return g * std::pow(std::numbers::pi, 1.5) / (4 * eta * eta)
                   * (eta * r_i * r_j * helpers::s3(t)
                      + ((i == j) + 2 * alpha * helpers::component_difference(a, b, j) * r_i -
                         2 * beta * helpers::component_difference(a, b, i) * r_j) * helpers::s2(t)
                      + (2 * (i == j) - 4 * alpha * beta * helpers::component_difference(a, b, i) *
                                        helpers::component_difference(a, b, j) / eta) * helpers::s1(t));
        });
    });
}

Eigen::MatrixXd overlap(const BasisSet& basis_set) {
    return one_electron_matrix(basis_set, OverlapIntegrals());
}

Eigen::MatrixXd kinetic_energy(const BasisSet& basis_set) {
    return one_electron_matrix(basis_set, KineticEnergyIntegrals());
}

Eigen::MatrixXd potential_energy(const BasisSet& basis_set, const Molecule& mol) {
    return one_electron_matrix(basis_set, NuclearPotentialIntegrals(mol));
}
