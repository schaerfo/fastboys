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

#include "DensityMatrix.hpp"

void DensityMatrix::updateDensity(const Eigen::MatrixXd& c) {
    previous_density_ = static_cast<Eigen::MatrixXd>(*this);
    for (Eigen::MatrixXd::Index mu=0; mu<rows(); ++mu){
        for (Eigen::MatrixXd::Index nu=0; nu<cols(); ++nu){
            double currValue = 0;
            for (decltype(n_occ_) m=0; m<n_occ_; ++m){
                currValue += c(mu, m) * c(nu, m);
            }
            (*this)(mu, nu) = 2 * currValue;
        }
    }
}

double DensityMatrix::difference_norm() const {
    return (*this - previous_density_).norm() / cols();
}

void DampingDensityMatrix::updateDensity(const Eigen::MatrixXd& c) {
    DensityMatrix::updateDensity(c);
    if (current_iteration < n_) {
        ++current_iteration;
        static_cast<Eigen::MatrixXd&>(*this) = (1 - alpha_) * *this + alpha_ * previous_density_;
    }
}
