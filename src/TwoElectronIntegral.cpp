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

#include "TwoElectronIntegral.hpp"
#include "TwoElectronIntegralDetail.hpp"

#include "Helpers.hpp"

#include <algorithm>

#ifdef ENABLE_OPENMP
# include <omp.h>
#endif

using namespace helpers;

void order_orbital_tuple(const BasisFunction*& r, const BasisFunction*& s, const BasisFunction*& t, const BasisFunction*& u) {
    if (is_s_orbital(r->type) && is_p_orbital(s->type))
        std::swap(r, s);

    if (is_s_orbital(t->type) && is_p_orbital(u->type))
        std::swap(t, u);

    if ((is_s_orbital(s->type) && is_p_orbital(t->type) && is_p_orbital(u->type)) ||
        (is_s_orbital(r->type) && is_s_orbital(s->type) && is_p_orbital(t->type))) {
        std::swap(r, t);
        std::swap(s, u);
    }
}

template <typename F>
double add_primitive_integrals(const BasisFunction& r, const BasisFunction& s, const BasisFunction& t, const BasisFunction& u, F primitive_func) {
    double res = 0;
    for (const auto& [alpha, a]: r.primitives) {
        for (const auto&[beta, b]: s.primitives) {
            for (const auto&[gamma, c]: t.primitives) {
                for (const auto&[delta, d]: u.primitives) {
                    res += a * b * c * d * primitive_func(alpha, beta, gamma, delta);
                }
            }
        }
    }
    return res;
}

double prim_ssss(double alpha, double beta, double gamma, double delta, const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d) {
    auto eta = alpha + beta, theta = gamma + delta;
    auto q = std::sqrt(eta*theta / (eta+theta));
    auto r = (alpha * a + beta * b) / eta - (gamma * c + delta * d) / theta;
    auto t = q * r.norm();
    return gaussian_product(alpha, a, beta, b) * gaussian_product(gamma, c, delta, d) * q * pi_cubed * s1(t) / std::pow(eta * theta, 3. / 2.);
}

double coulomb_ssss(const BasisFunction& r, const BasisFunction& s, const BasisFunction& t, const BasisFunction& u) {
    assert(is_s_orbital(r.type) && is_s_orbital(s.type) && is_s_orbital(t.type) && is_s_orbital(u.type));
    return add_primitive_integrals(r, s, t, u, [a = r.center, b = s.center, c = t.center, d = u.center](double alpha, double beta, double gamma, double delta) {
        return prim_ssss(alpha, beta, gamma, delta, a, b, c, d);
    });
}

double prim_psss(double alpha, double beta, double gamma, double delta, const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d, BasisFunction::Type i) {
    auto eta = alpha + beta, theta = gamma + delta;
    auto q = std::sqrt(eta*theta / (eta+theta));
    auto r = (alpha * a + beta * b) / eta - (gamma * c + delta * d) / theta;
    auto t = q * r.norm();
    return gaussian_product(alpha, a, beta, b) * gaussian_product(gamma, c, delta, d) * q * pi_cubed * (q*q * vector_component(r, i) * s2(t) - 2 * beta * component_difference(a, b, i) * s1(t))
                                                                                            / (2 * std::pow(eta, 5./2.) * std::pow(theta, 3./2.));
}

double coulomb_psss(const BasisFunction& r, const BasisFunction& s, const BasisFunction& t, const BasisFunction& u) {
    assert(is_p_orbital(r.type) && is_s_orbital(s.type) && is_s_orbital(t.type) && is_s_orbital(u.type));
    return add_primitive_integrals(r, s, t, u, [a = r.center, i = r.type, b = s.center, c = t.center, d = u.center](double alpha, double beta, double gamma, double delta) {
        return prim_psss(alpha, beta, gamma, delta, a, b, c, d, i);
    });
}

double prim_ppss(double alpha, double beta, double gamma, double delta, const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d, BasisFunction::Type i, BasisFunction::Type j) {
    auto delta_ab_i = component_difference(a, b, i);
    auto delta_ab_j = component_difference(a, b, j);

    auto eta = alpha + beta, theta = gamma + delta;
    auto q = std::sqrt(eta*theta / (eta+theta));
    auto r = (alpha * a + beta * b) / eta - (gamma * c + delta * d) / theta;
    auto r_i = vector_component(r, i);
    auto r_j = vector_component(r, j);
    auto t = q * r.norm();
    return gaussian_product(alpha, a, beta, b) * gaussian_product(gamma, c, delta, d) * q * pi_cubed/ (4 * std::pow(eta, 3.5) * std::pow(theta, 1.5)) *
            (q * q * q * q * r_i * r_j * s3(t)
             + q * q * ((i == j) + 2 * alpha * delta_ab_j * r_i - 2 * beta * delta_ab_i * r_j) * s2(t)
             + (2 * eta * (i == j) - 4 * alpha * beta * delta_ab_i * delta_ab_j) * s1(t));
}

double coulomb_ppss(const BasisFunction& r, const BasisFunction& s, const BasisFunction& t, const BasisFunction& u) {
    assert(is_p_orbital(r.type) && is_p_orbital(s.type) && is_s_orbital(t.type) && is_s_orbital(u.type));
    return add_primitive_integrals(r, s, t, u, [a = r.center, i = r.type, b = s.center, j = s.type, c = t.center, d = u.center](double alpha, double beta, double gamma, double delta) {
        return prim_ppss(alpha, beta, gamma, delta, a, b, c, d, i, j);
    });
}

double prim_psps(double alpha, double beta, double gamma, double delta, const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d, BasisFunction::Type i, BasisFunction::Type k) {
    auto delta_ab_i = component_difference(a, b, i);
    auto delta_cd_k = component_difference(c, d, k);

    auto eta = alpha + beta, theta = gamma + delta;
    auto q = std::sqrt(eta*theta / (eta+theta));
    auto r = (alpha * a + beta * b) / eta - (gamma * c + delta * d) / theta;
    auto r_i = vector_component(r, i);
    auto r_k = vector_component(r, k);
    auto t = q * r.norm();
    return gaussian_product(alpha, a, beta, b) * gaussian_product(gamma, c, delta, d) * q * pi_cubed/ (4 * std::pow(eta, 2.5) * std::pow(theta, 2.5)) *
           (-q * q * q * q * r_i * r_k * s3(t)
            + q * q * (2 * beta * delta_ab_i * r_k - 2 * delta * delta_cd_k * r_i - (i == k)) * s2(t)
            + 4 * beta * delta * delta_ab_i * delta_cd_k * s1(t));
}

double coulomb_psps(const BasisFunction& r, const BasisFunction& s, const BasisFunction& t, const BasisFunction& u) {
    assert(is_p_orbital(r.type) && is_s_orbital(s.type) && is_p_orbital(t.type) && is_s_orbital(u.type));
    return add_primitive_integrals(r, s, t, u, [a = r.center, i = r.type, b = s.center, c = t.center, k = t.type, d = u.center](double alpha, double beta, double gamma, double delta) {
        return prim_psps(alpha, beta, gamma, delta, a, b, c, d, i, k);
    });
}

double prim_ppps(double alpha, double beta, double gamma, double delta, const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d, BasisFunction::Type i, BasisFunction::Type j, BasisFunction::Type k) {
    auto delta_ij = i == j;
    auto delta_jk = j == k;

    auto delta_ab_i = component_difference(a, b, i);
    auto delta_ab_j = component_difference(a, b, j);
    auto delta_cd_k = component_difference(c, d, k);

    auto eta = alpha + beta, theta = gamma + delta;
    auto q = std::sqrt(eta*theta / (eta+theta));
    auto r = (alpha * a + beta * b) / eta - (gamma * c + delta * d) / theta;
    auto r_i = vector_component(r, i);
    auto r_j = vector_component(r, j);
    auto r_k = vector_component(r, k);
    auto t = q * r.norm();

    auto u1 = -4 * eta * delta * delta_cd_k * delta_ij;
    auto u2 = 2 * (2 * beta * delta * delta_ab_i * delta_cd_k * r_j + beta * delta_ab_i * delta_jk - eta * r_k * delta_ij - delta * delta_cd_k * delta_ij);
    auto u3 = 2 * beta * delta_ab_i * r_j * r_k - 2 * delta * delta_cd_k * r_i * r_j - r_k * delta_ij - r_i * delta_jk - r_j * (i == k);
    auto u4 = -r_i * r_j * r_k;

    return alpha * delta_ab_j * prim_psps(alpha, beta, gamma, delta, a, b, c, d, i, k) / eta
           + gaussian_product(alpha, a, beta, b) * gaussian_product(gamma, c, delta, d) * q * pi_cubed/ (8 * std::pow(eta, 3.5) * std::pow(theta, 2.5)) *
           (q * q * q * q * q * q * u4 * s4(t) + q * q * q * q * u3 * s3(t) + q * q * u2 * s2(t) + u1 * s1(t));
}

double coulomb_ppps(const BasisFunction& r, const BasisFunction& s, const BasisFunction& t, const BasisFunction& u) {
    assert(is_p_orbital(r.type) && is_p_orbital(s.type) && is_p_orbital(t.type) && is_s_orbital(u.type));
    return add_primitive_integrals(r, s, t, u, [a = r.center, i = r.type, b = s.center, j = s.type, c = t.center, k = t.type, d = u.center](double alpha, double beta, double gamma, double delta) {
        return prim_ppps(alpha, beta, gamma, delta, a, b, c, d, i, j, k);
    });
}

double prim_pppp(double alpha, double beta, double gamma, double delta, const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& d, BasisFunction::Type i, BasisFunction::Type j, BasisFunction::Type k, BasisFunction::Type l) {
    auto delta_ij = i == j;
    auto delta_ik = i == k;
    auto delta_il = i == l;
    auto delta_jk = j == k;
    auto delta_jl = j == l;
    auto delta_kl = k == l;

    auto delta_ab_i = component_difference(a, b, i);
    auto delta_ab_j = component_difference(a, b, j);
    auto delta_cd_k = component_difference(c, d, k);
    auto delta_cd_l = component_difference(c, d, l);

    auto eta = alpha + beta, theta = gamma + delta;
    auto q = std::sqrt(eta*theta / (eta+theta));
    auto r = (alpha * a + beta * b) / eta - (gamma * c + delta * d) / theta;
    auto r_i = vector_component(r, i);
    auto r_j = vector_component(r, j);
    auto r_k = vector_component(r, k);
    auto r_l = vector_component(r, l);
    auto t = q * r.norm();

    auto u1 = -4 * eta * delta * delta_cd_k * delta_ij;
    auto u2 = 2 * (2 * beta * delta * delta_ab_i * delta_cd_k * r_j + beta * delta_ab_i * delta_jk - eta * r_k * delta_ij - delta * delta_cd_k * delta_ij);
    auto u3 = 2 * beta * delta_ab_i * r_j * r_k - 2 * delta * delta_cd_k * r_i * r_j - r_k * delta_ij - r_i * delta_jk - r_j * (i == k);
    auto u4 = -r_i * r_j * r_k;

    auto v1 = 4 * eta * theta * delta_ij * delta_kl;
    auto v2 = 2 * ((eta + theta) * delta_ij * delta_kl - 2 * beta * theta * delta_ab_i * r_j * delta_kl - 2 * beta * delta * delta_ab_i * delta_cd_k * delta_jl);
    auto v3 = 2 * theta * r_i * r_j * delta_kl + 2 * delta * delta_cd_k * (r_i * delta_jl + r_j * delta_il)
              - 2 * beta * delta_ab_i * (r_k * delta_jl + r_j * delta_kl) + delta_ij * delta_kl + delta_il * delta_jk + delta_ik * delta_jl;
    auto v4 = r_i * r_k * delta_jl + r_i * r_j * delta_kl + r_j * r_k * delta_il;

    return alpha * delta_ab_j * prim_ppps(gamma, delta, alpha, beta, c, d, a, b, k, l, i) / eta
           + gaussian_product(alpha, a, beta, b) * gaussian_product(gamma, c, delta, d) * q * pi_cubed/ (16 * std::pow(eta, 3.5) * std::pow(theta, 3.5)) *
           (-q * q * q * q * q * q * q * q * u4 * r_l * s5(t)
            + q * q * q * q * q * q * (v4 + 2 * gamma * delta_cd_l * u4 - u3 * r_l) * s4(t)
            + q * q * q * q * (v3 + 2 * gamma * delta_cd_l * u3 - u2 * r_l) * s3(t)
            + q * q * (v2 + 2 * gamma * delta_cd_l * u2 - u1 * r_l) * s2(t)
            + (v1 + 2 * gamma * delta_cd_l * u1) * s1(t));
}

double coulomb_pppp(const BasisFunction& r, const BasisFunction& s, const BasisFunction& t, const BasisFunction& u) {
    assert(is_p_orbital(r.type) && is_p_orbital(s.type) && is_p_orbital(t.type) && is_p_orbital(u.type));
    return add_primitive_integrals(r, s, t, u, [a = r.center, i = r.type, b = s.center, j = s.type, c = t.center, k = t.type, d = u.center, l = u.type](double alpha, double beta, double gamma, double delta) {
        return prim_pppp(alpha, beta, gamma, delta, a, b, c, d, i, j, k, l);
    });
}

std::vector<TwoElectronIntegral> calculate_two_electron_integrals(const BasisSet& basis_set) {
    auto n = basis_set.size();

#ifdef ENABLE_OPENMP
    std::vector<std::vector<TwoElectronIntegral>> thread_results(n);
    #pragma omp parallel for schedule(dynamic, 1) default(none) shared(thread_results, n, basis_set)
    for (std::ptrdiff_t mu = n-1; mu >= 0; --mu) {
        auto& result = thread_results[mu];
#else
    std::vector<TwoElectronIntegral> result;
    for (uint16_t mu=0; mu<n; ++mu) {
#endif
        for (uint16_t  nu=0; nu<=mu; ++nu) {
            for (uint16_t  lambda=0; lambda<=mu; ++lambda) {
                for (uint16_t  sigma=0; sigma<=(lambda == mu ? nu : lambda); ++sigma) {
                    const auto* r = &basis_set[mu];
                    const auto* s = &basis_set[nu];
                    const auto* t = &basis_set[lambda];
                    const auto* u = &basis_set[sigma];

                    order_orbital_tuple(r, s, t, u);

                    auto curr_integral = [&r = *r, &s = *s, &t = *t, &u = *u]() {
                        if (is_p_orbital(u.type))
                            return coulomb_pppp(r, s, t, u);
                        else if (is_p_orbital(t.type) && is_p_orbital(s.type))
                            return coulomb_ppps(r, s, t, u);
                        else if (is_p_orbital(s.type))
                            return coulomb_ppss(r, s, t, u);
                        else if (is_p_orbital(t.type))
                            return coulomb_psps(r, s, t, u);
                        else if (is_p_orbital(r.type))
                            return coulomb_psss(r, s, t, u);
                        else
                            return coulomb_ssss(r, s, t, u);
                    }();

                    if (std::abs(curr_integral) > 1e-9)
                        result.emplace_back(TwoElectronIntegral{static_cast<uint16_t>(mu), nu, lambda, sigma, curr_integral});
                }
            }
        }
#ifdef ENABLE_OPENMP
        result.shrink_to_fit();
#endif
    }

#ifdef ENABLE_OPENMP
    auto view = thread_results | std::ranges::views::transform([](const auto& vec){return vec.size();});

    std::vector<TwoElectronIntegral> result(std::accumulate(view.begin(), view.end(), std::size_t(0)));
    auto it = result.begin();
    for (const auto& curr_vector: thread_results)
        it = std::copy(curr_vector.begin(), curr_vector.end(), it);
#endif
    result.shrink_to_fit();
    return result;
}
