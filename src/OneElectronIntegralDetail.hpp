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

#ifndef ONEELECTRONINTEGRALDETAIL_HPP
#define ONEELECTRONINTEGRALDETAIL_HPP

struct OverlapIntegrals {
    double integral_ss(const BasisFunction& r, const BasisFunction& s);
    double integral_ps(const BasisFunction& r, const BasisFunction& s);
    double integral_pp(const BasisFunction& r, const BasisFunction& s);
};

struct KineticEnergyIntegrals {
    double integral_ss(const BasisFunction& r, const BasisFunction& s);
    double integral_ps(const BasisFunction& r, const BasisFunction& s);
    double integral_pp(const BasisFunction& r, const BasisFunction& s);
};

struct NuclearPotentialIntegrals {
    NuclearPotentialIntegrals(const Molecule& mol): mol_(mol) {}

    double integral_ss(const BasisFunction& r, const BasisFunction& s);
    double integral_ps(const BasisFunction& r, const BasisFunction& s);
    double integral_pp(const BasisFunction& r, const BasisFunction& s);

private:
    const Molecule& mol_;

    template<typename Func>
    double nuclear_loop(Func f) {
        double res = 0;
        for (const auto& curr_atom: mol_.get_atoms()) {
            res -= curr_atom.charge_ * f(curr_atom.pos_);
        }
        return res;
    }
};

#endif //ONEELECTRONINTEGRALDETAIL_HPP
