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

#include <gtest/gtest.h>

#include <string>
#include <sstream>

#include "CommonTestingData.hpp"
#include "BasisSet.hpp"

TEST(MoleculeTest, read_xyz_file_test) {
    std::string water_file(" 3\n"
                           " water\n"
                           "O 0.00000    0.00000     0.0000 \n"
                           "H 0.75690    0.58585     0.0000\n"
                           "H -0.75690    0.58585     0.0000");

    std::istringstream water_stream(water_file);
    Molecule water_test{water_stream};
    EXPECT_EQ(water_test, water);

    std::string ethene_file("6\n"
                            "\n"
                            "C         -3.18633        1.82554        0.00291\n"
                            "C         -2.24348        0.87638        0.04828\n"
                            "H         -2.93436        2.87687       -0.09562\n"
                            "H         -4.23885        1.56591        0.06285\n"
                            "H         -1.24643        1.13630       -0.00967\n"
                            "H         -2.51365       -0.11528        0.14121");
    std::istringstream ethene_stream(ethene_file);
    Molecule ethene_test(ethene_stream);
    EXPECT_EQ(ethene_test, ethene);

    std::string bad_file("1\n\nFs 0 0 0");
    std::istringstream bad_stream(bad_file);
    EXPECT_THROW(Molecule{bad_stream}, std::runtime_error);
}

TEST(MoleculeTest, construct_basisset_test) {
    auto basisset = nlohmann::json::parse(basisset_file);
    auto water_test = water.construct_basis_set(basisset);
    EXPECT_EQ(water_test, water_basis);
    auto ethene_test = ethene.construct_basis_set(basisset);
    EXPECT_EQ(ethene_test, ethene_basis);

    Molecule m({Atom("Sc", {0, 0, 0})});
    EXPECT_THROW(m.construct_basis_set(basisset), std::runtime_error);
}

TEST(MoleculeTest, atom_equals_test) {
    Atom a("H", {0, 0, 0});
    Atom b = a;
    EXPECT_EQ(a, b);
    b.pos_[1] = 0.5;
    EXPECT_NE(a, b);
    b.charge_ = 2;
    EXPECT_NE(a, b);
    b.pos_ = a.pos_;
    EXPECT_NE(a, b);
}
