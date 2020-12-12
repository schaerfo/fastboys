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

#ifndef DENSITYMATRIX_HPP
#define DENSITYMATRIX_HPP

#include <memory>

#include <Eigen/Dense>

class DensityMatrix: public Eigen::MatrixXd {
public:
    virtual ~DensityMatrix() = default;

    DensityMatrix(unsigned n, unsigned n_occ):
        Eigen::MatrixXd(Eigen::MatrixXd::Zero(n, n)),
        previous_density_(Eigen::MatrixXd::Zero(n, n)),
        n_occ_(n_occ)
    {}

    virtual void updateDensity(const Eigen::MatrixXd& c);
    double difference_norm() const;

protected:
    Eigen::MatrixXd previous_density_;

private:
    unsigned n_occ_;
};

class DampingDensityMatrix: public DensityMatrix {
public:
    DampingDensityMatrix(double alpha, unsigned n_damp, unsigned n, unsigned n_occ):
        DensityMatrix(n, n_occ),
        alpha_(alpha),
        n_(n_damp)
    {}

    void updateDensity(const Eigen::MatrixXd& c) override;

private:
    double alpha_;
    unsigned n_, current_iteration = 0;
};

class DensityMatrixFactory {
public:
    virtual ~DensityMatrixFactory() = default;

    virtual std::unique_ptr<DensityMatrix> get_density_matrix(unsigned n, unsigned n_occ) const {
        return std::make_unique<DensityMatrix>(n, n_occ);
    }
};

class DampingDensityMatrixFactory: public DensityMatrixFactory {
public:
    DampingDensityMatrixFactory(double alpha, unsigned int n_damp):
        alpha_(alpha),
        n_damp_(n_damp)
    {}

    std::unique_ptr<DensityMatrix> get_density_matrix(unsigned int n, unsigned int n_occ) const override {
        return std::make_unique<DampingDensityMatrix>(alpha_, n_damp_, n, n_occ);
    }

private:
    double alpha_;
    unsigned n_damp_;
};



#endif //DENSITYMATRIX_HPP
