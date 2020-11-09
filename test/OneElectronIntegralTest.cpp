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

#include "Molecule.hpp"
#include "OneElectronIntegrals.hpp"
#include "OneElectronIntegralDetail.hpp"

class OneElectronIntegralTest: public ::testing::Test {
protected:
    void SetUp() override {
        using PrimFunc = BasisFunction::PrimitiveFunction;

        Eigen::Vector3d water_o_pos = {0, 0, 0};
        Eigen::Vector3d water_h1_pos = {1.430541, 1.1072565, 0};
        Eigen::Vector3d water_h2_pos = {-1.430541, 1.1072565, 0};

        water.setAtoms({{Atom("O", water_o_pos)}, Atom{"H", water_h1_pos}, Atom{"H", water_h2_pos}});
        water_basis = {
                BasisFunction{water_o_pos, {PrimFunc{130.7093214, 4.251943277890349}, PrimFunc{23.80886605, 4.112294424303829}, PrimFunc{6.443608313, 1.2816225514654456}}, BasisFunction::Type::s},
                BasisFunction{water_o_pos, {PrimFunc{5.033151319, -0.2394130049495816}, PrimFunc{1.169596125, 0.3202342354447327}, PrimFunc{0.38038896, 0.24168555456024995}}, BasisFunction::Type::s},
                BasisFunction{water_o_pos, {PrimFunc{5.033151319, 1.6754501961028012}, PrimFunc{1.169596125, 1.0535680440009998}, PrimFunc{0.38038896, 0.1669028790864184}}, BasisFunction::Type::px},
                BasisFunction{water_o_pos, {PrimFunc{5.033151319, 1.6754501961028012}, PrimFunc{1.169596125, 1.0535680440009998}, PrimFunc{0.38038896, 0.1669028790864184}}, BasisFunction::Type::py},
                BasisFunction{water_o_pos, {PrimFunc{5.033151319, 1.6754501961028012}, PrimFunc{1.169596125, 1.0535680440009998}, PrimFunc{0.38038896, 0.1669028790864184}}, BasisFunction::Type::pz},
                BasisFunction{water_h1_pos, {PrimFunc{3.425250914, 0.2769343550887434}, PrimFunc{0.6239137298, 0.26783885161885207}, PrimFunc{0.168855404, 0.08347367113276243}}, BasisFunction::Type::s},
                BasisFunction{water_h2_pos, {PrimFunc{3.425250914, 0.2769343550887434}, PrimFunc{0.6239137298, 0.26783885161885207}, PrimFunc{0.168855404, 0.08347367113276243}}, BasisFunction::Type::s}
        };

        Eigen::Vector3d ethene_c1_pos = {-6.0221637, 3.4502706, 0.0054999};
        Eigen::Vector3d ethene_c2_pos = {-4.2401772, 1.6563582, 0.0912492};
        Eigen::Vector3d ethene_h1_pos = {-5.5459404, 5.4372843, -0.1807218};
        Eigen::Vector3d ethene_h2_pos = {-8.0114265, 2.9595699, 0.1187865};
        Eigen::Vector3d ethene_h3_pos = {-2.3557527, 2.147607, -0.0182763};
        Eigen::Vector3d ethene_h4_pos = {-4.7507985, -0.2178792, 0.2668869};

        ethene.setAtoms({{Atom("C", ethene_c1_pos)}, {Atom("C", ethene_c2_pos)}, {Atom("H", ethene_h1_pos)}, {Atom("H", ethene_h2_pos)}, {Atom("H", ethene_h3_pos)}, {Atom("H", ethene_h4_pos)}});

        ethene_basis = {
                BasisFunction{ethene_c1_pos, {PrimFunc{71.61683735, 2.707814416331934}, PrimFunc{13.04509632, 2.6188802153557194}, PrimFunc{3.53051216, 0.8161905733612459}}, BasisFunction::Type::s},
                BasisFunction{ethene_c1_pos, {PrimFunc{2.941249355, -0.1600171867314577}, PrimFunc{0.6834830964, 0.21403591446739756}, PrimFunc{0.2222899159, 0.16153609753464734}}, BasisFunction::Type::s},
                BasisFunction{ethene_c1_pos, {PrimFunc{2.941249355, 0.8560445064070002}, PrimFunc{0.6834830964, 0.5383037575547593}, PrimFunc{0.2222899159, 0.08527635921608212}}, BasisFunction::Type::px},
                BasisFunction{ethene_c1_pos, {PrimFunc{2.941249355, 0.8560445064070002}, PrimFunc{0.6834830964, 0.5383037575547593}, PrimFunc{0.2222899159, 0.08527635921608212}}, BasisFunction::Type::py},
                BasisFunction{ethene_c1_pos, {PrimFunc{2.941249355, 0.8560445064070002}, PrimFunc{0.6834830964, 0.5383037575547593}, PrimFunc{0.2222899159, 0.08527635921608212}}, BasisFunction::Type::pz},
                BasisFunction{ethene_c2_pos, {PrimFunc{71.61683735, 2.707814416331934}, PrimFunc{13.04509632, 2.6188802153557194}, PrimFunc{3.53051216, 0.8161905733612459}}, BasisFunction::Type::s},
                BasisFunction{ethene_c2_pos, {PrimFunc{2.941249355, -0.1600171867314577}, PrimFunc{0.6834830964, 0.21403591446739756}, PrimFunc{0.2222899159, 0.16153609753464734}}, BasisFunction::Type::s},
                BasisFunction{ethene_c2_pos, {PrimFunc{2.941249355, 0.8560445064070002}, PrimFunc{0.6834830964, 0.5383037575547593}, PrimFunc{0.2222899159, 0.08527635921608212}}, BasisFunction::Type::px},
                BasisFunction{ethene_c2_pos, {PrimFunc{2.941249355, 0.8560445064070002}, PrimFunc{0.6834830964, 0.5383037575547593}, PrimFunc{0.2222899159, 0.08527635921608212}}, BasisFunction::Type::py},
                BasisFunction{ethene_c2_pos, {PrimFunc{2.941249355, 0.8560445064070002}, PrimFunc{0.6834830964, 0.5383037575547593}, PrimFunc{0.2222899159, 0.08527635921608212}}, BasisFunction::Type::pz},
                BasisFunction{ethene_h1_pos, {PrimFunc{3.425250914, 0.2769343550887434}, PrimFunc{0.6239137298, 0.26783885161885207}, PrimFunc{0.168855404, 0.08347367113276243}}, BasisFunction::Type::s},
                BasisFunction{ethene_h2_pos, {PrimFunc{3.425250914, 0.2769343550887434}, PrimFunc{0.6239137298, 0.26783885161885207}, PrimFunc{0.168855404, 0.08347367113276243}}, BasisFunction::Type::s},
                BasisFunction{ethene_h3_pos, {PrimFunc{3.425250914, 0.2769343550887434}, PrimFunc{0.6239137298, 0.26783885161885207}, PrimFunc{0.168855404, 0.08347367113276243}}, BasisFunction::Type::s},
                BasisFunction{ethene_h4_pos, {PrimFunc{3.425250914, 0.2769343550887434}, PrimFunc{0.6239137298, 0.26783885161885207}, PrimFunc{0.168855404, 0.08347367113276243}}, BasisFunction::Type::s}
        };
    }

    Molecule water, ethene;
    BasisSet water_basis, ethene_basis;
};

TEST_F(OneElectronIntegralTest, overlap_matrix_test) {
    Eigen::MatrixXd water_expected(water_basis.size(), water_basis.size());
    water_expected <<     1,  0.236704,         0,         0,         0, 0.0539675, 0.0539675,
                   0.236704,         1,         0,         0,         0,  0.474755,  0.474755,
                          0,         0,         1,         0,         0,  0.311117, -0.311117,
                          0,         0,         0,         1,         0,  0.240809,  0.240809,
                          0,         0,         0,         0,         1,         0,         0,
                  0.0539675,  0.474755,  0.311117,  0.240809,         0,         1,  0.251674,
                  0.0539675,  0.474755, -0.311117,  0.240809,         0,  0.251674,         1;

    auto water_actual = overlap(water_basis);
    EXPECT_TRUE(water_expected.isApprox(water_actual, 1e-6));

    Eigen::MatrixXd ethene_expected(ethene_basis.size(), ethene_basis.size());
    ethene_expected <<
            1,     0.248362,            0,            0,            0,  2.45447e-06,    0.0433861,    -0.050683,    0.0510222,  -0.00243887,    0.0628978,    0.0628771,   0.00689784,    0.0068979,
            0.248362,            1,            0,            0,            0,    0.0433861,     0.395982,    -0.288586,     0.290517,   -0.0138868,     0.493207,     0.493126,     0.115778,     0.115779,
            0,            0,            1,            0,            0,     0.050683,     0.288586,   -0.0429629,     0.282607,   -0.0135086,     0.108796,    -0.454351,     0.137535,     0.047692,
            0,            0,            0,            1,            0,   -0.0510222,    -0.290517,     0.282607,    -0.046733,    0.0135991,     0.453943,    -0.112077,   -0.0488658,    -0.137601,
            0,            0,            0,            0,            1,   0.00243887,    0.0138868,   -0.0135086,    0.0135991,     0.237115,   -0.0425432,    0.0258748, -0.000891897,   0.00980527,
            2.45447e-06,    0.0433861,     0.050683,   -0.0510222,   0.00243887,            1,     0.248362,            0,            0,            0,   0.00590906,   0.00606029,    0.0712792,    0.0712798,
            0.0433861,     0.395982,     0.288586,    -0.290517,    0.0138868,     0.248362,            1,            0,            0,            0,     0.103512,     0.105419,     0.524501,     0.524503,
            -0.050683,    -0.288586,   -0.0429629,     0.282607,   -0.0135086,            0,            0,            1,            0,            0,   -0.0431416,    -0.127196,     0.469828,     -0.12731,
            0.0510222,     0.290517,     0.282607,    -0.046733,    0.0135991,            0,            0,            0,            1,            0,     0.124919,    0.0439544,     0.122479,    -0.467291,
            -0.00243887,   -0.0138868,   -0.0135086,    0.0135991,     0.237115,            0,            0,            0,            0,            1,  -0.00898575,  0.000928772,   -0.0273071,    0.0437905,
            0.0628978,     0.493207,     0.108796,     0.453943,   -0.0425432,   0.00590906,     0.103512,   -0.0431416,     0.124919,  -0.00898575,            1,     0.150204,    0.0567686,    0.0170912,
            0.0628771,     0.493126,    -0.454351,    -0.112077,    0.0258748,   0.00606029,     0.105419,    -0.127196,    0.0439544,  0.000928772,     0.150204,            1,    0.0173428,    0.0584586,
            0.00689784,     0.115778,     0.137535,   -0.0488658, -0.000891897,    0.0712792,     0.524501,     0.469828,     0.122479,   -0.0273071,    0.0567686,    0.0173428,            1,     0.167245,
            0.0068979,     0.115779,     0.047692,    -0.137601,   0.00980527,    0.0712798,     0.524503,     -0.12731,    -0.467291,    0.0437905,    0.0170912,    0.0584586,     0.167245,            1;

    auto ethene_actual = overlap(ethene_basis);
    EXPECT_TRUE(ethene_expected.isApprox(ethene_actual, 1e-6));
}

TEST_F(OneElectronIntegralTest, kinetic_energy_matrix_test) {
    Eigen::MatrixXd water_expected(water_basis.size(), water_basis.size());
    water_expected <<     29.0032,   -0.168011,           0,           0,           0, -0.00251709, -0.00251709,
                        -0.168011,    0.808128,           0,           0,           0,    0.128572,    0.128572,
                                0,           0,     2.52873,           0,           0,    0.224756,   -0.224756,
                                0,           0,           0,     2.52873,           0,    0.173964,    0.173964,
                                0,           0,           0,           0,     2.52873,           0,           0,
                      -0.00251709,    0.128572,    0.224756,    0.173964,           0,    0.760032,  0.00847527,
                      -0.00251709,    0.128572,   -0.224756,    0.173964,           0,  0.00847527,    0.760032;

    auto water_actual = kinetic_energy(water_basis);
    EXPECT_TRUE(water_expected.isApprox(water_actual, 1e-6));

    Eigen::MatrixXd ethene_expected(ethene_basis.size(), ethene_basis.size());
    ethene_expected <<
            15.8911,     -0.08589,            0,            0,            0, -8.52904e-05,  -0.00831933,   0.00701355,  -0.00706049,  0.000337493,   -0.0103719,   -0.0103771,  -0.00257274,  -0.00257275,
            -0.08589,      0.47225,            0,            0,            0,  -0.00831933,    0.0581925,    -0.120767,     0.121575,  -0.00581131,     0.108257,      0.10821,   -0.0141934,   -0.0141934,
            0,            0,      1.47773,            0,            0,  -0.00701355,     0.120767,    -0.125113,     0.218422,   -0.0104406,      0.06123,     -0.25565,   0.00381032,   0.00132133,
            0,            0,            0,      1.47773,            0,   0.00706049,    -0.121575,     0.218422,    -0.128026,    0.0105105,     0.255478,   -0.0630623,   -0.0013538,  -0.00381232,
            0,            0,            0,            0,      1.47773, -0.000337493,   0.00581131,   -0.0104406,    0.0105105,     0.091355,   -0.0239433,     0.014559, -2.47095e-05,   0.00027166,
            -8.52904e-05,  -0.00831933,  -0.00701355,   0.00706049, -0.000337493,      15.8911,     -0.08589,            0,            0,            0,  -0.00241289,  -0.00243803,  -0.00774103,   -0.0077408,
            -0.00831933,    0.0581925,     0.120767,    -0.121575,   0.00581131,     -0.08589,      0.47225,            0,            0,            0,   -0.0147621,   -0.0146924,     0.126753,     0.126755,
            0.00701355,    -0.120767,    -0.125113,     0.218422,   -0.0104406,            0,            0,      1.47773,            0,            0, -0.000294993,  -0.00129719,     0.287804,   -0.0779871,
            -0.00706049,     0.121575,     0.218422,    -0.128026,    0.0105105,            0,            0,            0,      1.47773,            0,  0.000854173,  0.000448264,    0.0750274,    -0.286252,
            0.000337493,  -0.00581131,   -0.0104406,    0.0105105,     0.091355,            0,            0,            0,            0,      1.47773, -6.14427e-05,  9.47197e-06,   -0.0167276,    0.0268251,
            -0.0103719,     0.108257,      0.06123,     0.255478,   -0.0239433,  -0.00241289,   -0.0147621, -0.000294993,  0.000854173, -6.14427e-05,     0.760032,  -0.00764194,   -0.0110525,  -0.00637142,
            -0.0103771,      0.10821,     -0.25565,   -0.0630623,     0.014559,  -0.00243803,   -0.0146924,  -0.00129719,  0.000448264,  9.47197e-06,  -0.00764194,     0.760032,  -0.00642915,   -0.0111231,
            -0.00257274,   -0.0141934,   0.00381032,   -0.0013538, -2.47095e-05,  -0.00774103,     0.126753,     0.287804,    0.0750274,   -0.0167276,   -0.0110525,  -0.00642915,     0.760032,   -0.0058082,
            -0.00257275,   -0.0141934,   0.00132133,  -0.00381232,   0.00027166,   -0.0077408,     0.126755,   -0.0779871,    -0.286252,    0.0268251,  -0.00637142,   -0.0111231,   -0.0058082,     0.760032;

    auto ethene_actual = kinetic_energy(ethene_basis);
    std::cout << ethene_actual - ethene_expected << '\n';
    EXPECT_TRUE(ethene_expected.isApprox(ethene_actual, 1e-5));
}

TEST_F(OneElectronIntegralTest, potential_energy_matrix_test) {
    Eigen::MatrixXd water_expected(water_basis.size(), water_basis.size());
    water_expected <<     -61.724,   -7.44478,          0, -0.0190003,          0,   -1.74672,   -1.74672,
                         -7.44478,   -10.1429,          0,  -0.223497,          0,   -3.86944,   -3.86944,
                                0,          0,   -10.1426,          0,          0,   -2.25477,    2.25477,
                       -0.0190003,  -0.223497,          0,     -10.08,          0,   -1.81819,   -1.81819,
                                0,          0,          0,          0,   -9.98633,          0,          0,
                         -1.74672,   -3.86944,   -2.25477,   -1.81819,          0,   -5.83641,    -1.6165,
                         -1.74672,   -3.86944,    2.25477,   -1.81819,          0,    -1.6165,   -5.83641;

    auto water_actual = potential_energy(water_basis, water);
    EXPECT_TRUE(water_expected.isApprox(water_actual, 1e-5));

    Eigen::MatrixXd ethene_expected(ethene_basis.size(), ethene_basis.size());
    ethene_expected <<
            -37.513,     -5.15097,   -0.0408824,    0.0413713,  -0.00196553, -2.67482e-05,    -0.918897,      1.07502,     -1.08221,    0.0517299,     -1.32226,     -1.32193,    -0.145001,    -0.145005,
            -5.15097,     -9.00149,    -0.468939,     0.474287,    -0.022549,    -0.920344,     -3.61042,      2.73024,     -2.74793,     0.131355,     -3.73713,     -3.74031,    -0.948304,    -0.948507,
            -0.0408824,    -0.468939,     -9.08815,     0.202988,  -0.00259039,      -1.0767,     -2.73375,     0.707631,     -2.61236,     0.126028,    -0.897004,      3.09217,     -1.08772,    -0.465963,
            0.0413713,     0.474287,     0.202988,     -9.09179,    0.0263492,      1.08391,       2.7528,     -2.61335,     0.743021,    -0.123085,     -3.08624,     0.923586,     0.475139,      1.08944,
            -0.00196553,    -0.022549,  -0.00259039,    0.0263492,     -8.69349,   -0.0518112,    -0.131577,     0.126068,    -0.123078,      -1.8502,     0.295146,    -0.189528, -0.000369271,   -0.0744214,
            -2.67482e-05,    -0.920344,      -1.0767,      1.08391,   -0.0518112,     -37.5497,     -5.16009,    0.0398586,   -0.0400551,   0.00191471,    -0.124247,    -0.127454,     -1.50495,     -1.50497,
            -0.918897,     -3.61042,     -2.73375,       2.7528,    -0.131577,     -5.16009,     -9.02887,     0.462191,    -0.464467,    0.0222026,    -0.840759,    -0.855662,     -4.04999,     -4.05043,
            1.07502,      2.73024,     0.707631,     -2.61335,     0.126068,    0.0398586,     0.462191,     -9.12055,     0.190557,   -0.0011641,      0.41583,     0.993556,     -3.25734,      1.06011,
            -1.08221,     -2.74793,     -2.61236,     0.743021,    -0.123078,   -0.0400551,    -0.464467,     0.190557,     -9.12135,      0.02741,    -0.977157,    -0.423344,     -1.02564,      3.23857,
            0.0517299,     0.131355,     0.126028,    -0.123085,      -1.8502,   0.00191471,    0.0222026,   -0.0011641,      0.02741,     -8.71434,    0.0673712, -0.000543266,     0.204128,    -0.309931,
            -1.32226,     -3.73713,    -0.897004,     -3.08624,     0.295146,    -0.124247,    -0.840759,      0.41583,    -0.977157,    0.0673712,     -6.26052,     -1.03432,    -0.378049,    -0.139299,
            -1.32193,     -3.74031,      3.09217,     0.923586,    -0.189528,    -0.127454,    -0.855662,     0.993556,    -0.423344, -0.000543266,     -1.03432,     -6.26917,    -0.141347,    -0.388866,
            -0.145001,    -0.948304,     -1.08772,     0.475139, -0.000369271,     -1.50495,     -4.04999,     -3.25734,     -1.02564,     0.204128,    -0.378049,    -0.141347,     -6.44978,     -1.18658,
            -0.145005,    -0.948507,    -0.465963,      1.08944,   -0.0744214,     -1.50497,     -4.05043,      1.06011,      3.23857,    -0.309931,    -0.139299,    -0.388866,     -1.18658,     -6.45084;

    auto ethene_actual = potential_energy(ethene_basis, ethene);
    EXPECT_TRUE(ethene_expected.isApprox(ethene_actual, 1e-6));
}

class OneElectronIntegralDeathTest: public ::testing::Test {
protected:
    void SetUp() override {
        Eigen::Vector3d pos = {0, 0, 0};

        s = p = BasisFunction{pos, {}, BasisFunction::Type::s};
        p.type = BasisFunction::Type::px;
    }

    BasisFunction s, p;
    Molecule m;
};

TEST_F(OneElectronIntegralDeathTest, overlap_test) {
    OverlapIntegrals overlap;
    EXPECT_DEBUG_DEATH(overlap.integral_ss(s, p), "");
    EXPECT_DEBUG_DEATH(overlap.integral_ss(p, s), "");
    EXPECT_DEBUG_DEATH(overlap.integral_ss(p, p), "");
    EXPECT_DEBUG_DEATH(overlap.integral_ps(s, s), "");
    EXPECT_DEBUG_DEATH(overlap.integral_ps(s, p), "");
    EXPECT_DEBUG_DEATH(overlap.integral_ps(p, p), "");
    EXPECT_DEBUG_DEATH(overlap.integral_pp(s, s), "");
    EXPECT_DEBUG_DEATH(overlap.integral_pp(s, p), "");
    EXPECT_DEBUG_DEATH(overlap.integral_pp(p, s), "");
}

TEST_F(OneElectronIntegralDeathTest, kinetic_energy_test) {
    KineticEnergyIntegrals kinetic_energy;
    EXPECT_DEBUG_DEATH(kinetic_energy.integral_ss(s, p), "");
    EXPECT_DEBUG_DEATH(kinetic_energy.integral_ss(p, s), "");
    EXPECT_DEBUG_DEATH(kinetic_energy.integral_ss(p, p), "");
    EXPECT_DEBUG_DEATH(kinetic_energy.integral_ps(s, s), "");
    EXPECT_DEBUG_DEATH(kinetic_energy.integral_ps(s, p), "");
    EXPECT_DEBUG_DEATH(kinetic_energy.integral_ps(p, p), "");
    EXPECT_DEBUG_DEATH(kinetic_energy.integral_pp(s, s), "");
    EXPECT_DEBUG_DEATH(kinetic_energy.integral_pp(s, p), "");
    EXPECT_DEBUG_DEATH(kinetic_energy.integral_pp(p, s), "");
}

TEST_F(OneElectronIntegralDeathTest, nuclear_energy_test) {
    NuclearPotentialIntegrals nuclear_energy(m);
    EXPECT_DEBUG_DEATH(nuclear_energy.integral_ss(s, p), "");
    EXPECT_DEBUG_DEATH(nuclear_energy.integral_ss(p, s), "");
    EXPECT_DEBUG_DEATH(nuclear_energy.integral_ss(p, p), "");
    EXPECT_DEBUG_DEATH(nuclear_energy.integral_ps(s, s), "");
    EXPECT_DEBUG_DEATH(nuclear_energy.integral_ps(s, p), "");
    EXPECT_DEBUG_DEATH(nuclear_energy.integral_ps(p, p), "");
    EXPECT_DEBUG_DEATH(nuclear_energy.integral_pp(s, s), "");
    EXPECT_DEBUG_DEATH(nuclear_energy.integral_pp(s, p), "");
    EXPECT_DEBUG_DEATH(nuclear_energy.integral_pp(p, s), "");
}
