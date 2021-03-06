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

#include "FockMatrix.hpp"

#include <future>
#include <iostream>
#include <span>

double get_a(uint16_t i, uint16_t j, uint16_t k, uint16_t l) {
    if ((k == i || k == j) && l == k) return 0.5;
    else if (j != i && (k == i || k == j || l == i || l == j) && l != k) return 1.5;
    else if (k != i && k != j && l != i && l != j && l != k)  return 2;
    else return 1;
}

double get_b(uint16_t i, uint16_t j, uint16_t k, uint16_t l) {
    if (k == i && k != j && l != i && l != j) return -1;
    else if ((j == i || k == j ) && l == j) return 0.5;
    else if (k == i && (k == j || l == i) && l != j) return 1;
    else if (k != i && (k==j || l == i ) && l != j) return 1.5;
    else return -0.5;
}

Eigen::MatrixXd electron_repulsion_matrix_worker(std::span<const TwoElectronIntegral> span, const Eigen::MatrixXd& density) {
    auto size = density.cols();
    Eigen::MatrixXd result(size, size);
    result.setZero();
    for(auto [mu, nu, lambda, sigma, integral]: span) {
        // 1.
        result(mu, nu) += get_a(mu, nu, lambda, sigma) * density(lambda, sigma) * integral;

        // 2.
        if (lambda != nu)
            result(mu, lambda) += get_b(mu, nu, lambda, sigma) * density(nu, sigma) * integral;

        // 3.
        if (sigma != nu && sigma != lambda)
            result(mu, sigma) += get_b(mu, nu, sigma, lambda) * density(nu, lambda) * integral;

        // 4.
        if (nu >= lambda && nu != mu)
            result(nu, lambda) += get_b(nu, mu, lambda, sigma) * density(mu, sigma) * integral;

        // 5.
        if (nu < lambda && lambda != mu)
            result(lambda, nu) += get_b(lambda, sigma, nu, mu) * density(mu, sigma) * integral;

        // 6.
        if (nu >= sigma && sigma != lambda && nu != mu)
            result(nu, sigma) += get_b(nu, mu, sigma, lambda) * density(mu, lambda) * integral;

        // 7.
        if (nu < sigma && sigma != lambda)
            result(sigma, nu) += get_b(sigma, lambda, nu, mu) * density(mu, lambda) * integral;

        // 8.
        if (lambda != mu && lambda != nu && sigma != nu)
            result(lambda, sigma) += get_a(lambda, sigma, mu, nu) * density(mu, nu) * integral;
    }

    return result;
}

Eigen::MatrixXd electron_repulsion_matrix(const std::vector<TwoElectronIntegral>& integrals, const Eigen::MatrixXd& density, unsigned thread_count) {
    std::vector<std::future<Eigen::MatrixXd>> thread_results;
    auto n = integrals.size();

    for (unsigned i=0; i<thread_count; ++i) {
        auto start = static_cast<std::size_t>(i * static_cast<double>(n) / thread_count);
        auto end = static_cast<std::size_t>((i+1) * static_cast<double>(n) / thread_count);

        thread_results.emplace_back(std::async(electron_repulsion_matrix_worker, std::span(integrals.begin() + start, end - start), density));
    }

    auto size = density.cols();
    Eigen::MatrixXd result(size, size);
    result.setZero();

    for (auto& curr_result: thread_results)
        result += curr_result.get();

    for (Eigen::Index i=0; i<size; ++i) {
        for (Eigen::Index j=0; j<i; ++j) {
            result(j, i) = result(i, j);
        }
    }
    return result;
}

Eigen::MatrixXd initial_coefficients(const Eigen::MatrixXd& overlap) {
    Eigen::SelfAdjointEigenSolver<std::remove_reference_t<decltype(overlap)>> solver(overlap);
    const auto& eigvals = solver.eigenvalues();
    const auto& eigvecs = solver.eigenvectors();

    Eigen::MatrixXd result(overlap.rows(), overlap.cols());
    for (Eigen::MatrixXd::Index r=0; r<overlap.rows(); ++r){
        for (Eigen::MatrixXd::Index i=0; i<overlap.cols(); ++i){
            result(r, i) = eigvecs(r, i) / std::sqrt(eigvals(i));
        }
    }
    return result;
}

Eigen::MatrixXd update_coefficients(const Eigen::MatrixXd& c, const Eigen::MatrixXd& f_mo) {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(f_mo);
    return c * solver.eigenvectors();
}
