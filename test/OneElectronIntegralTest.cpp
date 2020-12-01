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

#include "IntegralTestData.hpp"

#include "CommonTestingData.hpp"
#include "OneElectronIntegrals.hpp"
#include "OneElectronIntegralDetail.hpp"


TEST(OneElectronIntegralTest, overlap_matrix_test) {
    Eigen::MatrixXd water_expected(water_basis.size(), water_basis.size());
    water_expected <<             1,    0.23670392,             0,             0,             0,    0.05396749,    0.05396749,
                         0.23670392,             1,             0,             0,             0,    0.47475471,    0.47475471,
                                  0,             0,             1,             0,             0,    0.31111725,   -0.31111725,
                                  0,             0,             0,             1,             0,    0.24080862,    0.24080862,
                                  0,             0,             0,             0,             1,             0,             0,
                         0.05396749,    0.47475471,    0.31111725,    0.24080862,             0,             1,    0.25167382,
                         0.05396749,    0.47475471,   -0.31111725,    0.24080862,             0,    0.25167382,             1;

    auto water_actual = overlap(water_basis);
    EXPECT_TRUE(water_expected.isApprox(water_actual, 1e-8));

    Eigen::MatrixXd ethene_expected(ethene_basis.size(), ethene_basis.size());
    ethene_expected <<
                    1,     0.2483624,             0,             0,             0, 2.4544732e-06,   0.043386084,  -0.050683007,   0.051022202, -0.0024388694,   0.062897846,   0.062877076,  0.0068978393,  0.0068978987,
            0.2483624,             1,             0,             0,             0,   0.043386084,     0.3959822,   -0.28858583,    0.29051719,  -0.013886768,    0.49320745,    0.49312589,    0.11577846,    0.11577919,
                    0,             0,             1,             0,             0,   0.050683007,    0.28858583,    -0.0429629,    0.28260668,  -0.013508645,    0.10879551,   -0.45435068,    0.13753509,   0.047692026,
                    0,             0,             0,             1,             0,  -0.051022202,   -0.29051719,    0.28260668,  -0.046733003,   0.013599051,    0.45394287,   -0.11207679,  -0.048865759,   -0.13760129,
                    0,             0,             0,             0,             1,  0.0024388694,   0.013886768,  -0.013508645,   0.013599051,    0.23711498,  -0.042543246,   0.025874833, -0.00089189723,  0.0098052673,
        2.4544732e-06,   0.043386084,   0.050683007,  -0.051022202,  0.0024388694,             1,     0.2483624,             0,             0,             0,  0.0059090599,  0.0060602869,   0.071279175,   0.071279815,
          0.043386084,     0.3959822,    0.28858583,   -0.29051719,   0.013886768,     0.2483624,             1,             0,             0,             0,    0.10351249,    0.10541882,    0.52450084,    0.52450312,
         -0.050683007,   -0.28858583,    -0.0429629,    0.28260668,  -0.013508645,             0,             0,             1,             0,             0,  -0.043141598,   -0.12719586,    0.46982774,   -0.12730971,
          0.051022202,    0.29051719,    0.28260668,  -0.046733003,   0.013599051,             0,             0,             0,             1,             0,    0.12491943,   0.043954435,    0.12247894,   -0.46729077,
        -0.0024388694,  -0.013886768,  -0.013508645,   0.013599051,    0.23711498,             0,             0,             0,             0,             1, -0.0089857515, 0.00092877194,  -0.027307073,   0.043790545,
          0.062897846,    0.49320745,    0.10879551,    0.45394287,  -0.042543246,  0.0059090599,    0.10351249,  -0.043141598,    0.12491943, -0.0089857515,             1,     0.1502037,   0.056768614,   0.017091161,
          0.062877076,    0.49312589,   -0.45435068,   -0.11207679,   0.025874833,  0.0060602869,    0.10541882,   -0.12719586,   0.043954435, 0.00092877194,     0.1502037,             1,   0.017342768,   0.058458587,
         0.0068978393,    0.11577846,    0.13753509,  -0.048865759, -0.00089189723,   0.071279175,    0.52450084,    0.46982774,    0.12247894,  -0.027307073,   0.056768614,   0.017342768,             1,    0.16724535,
         0.0068978987,    0.11577919,   0.047692026,   -0.13760129,  0.0098052673,   0.071279815,    0.52450312,   -0.12730971,   -0.46729077,   0.043790545,   0.017091161,   0.058458587,    0.16724535,             1;

    auto ethene_actual = overlap(ethene_basis);
    EXPECT_TRUE(ethene_expected.isApprox(ethene_actual, 1e-8));
}

TEST(OneElectronIntegralTest, kinetic_energy_matrix_test) {
    Eigen::MatrixXd water_expected(water_basis.size(), water_basis.size());
    water_expected <<     29.003204,   -0.16801096,             0,             0,             0,  -0.002517089,  -0.002517089,
                        -0.16801096,     0.8081279,             0,             0,             0,    0.12857248,    0.12857248,
                                  0,             0,     2.5287312,             0,             0,    0.22475609,   -0.22475609,
                                  0,             0,             0,     2.5287312,             0,    0.17396401,    0.17396401,
                                  0,             0,             0,             0,     2.5287312,             0,             0,
                       -0.002517089,    0.12857248,    0.22475609,    0.17396401,             0,    0.76003188,  0.0084752663,
                       -0.002517089,    0.12857248,   -0.22475609,    0.17396401,             0,  0.0084752663,    0.76003188;

    auto water_actual = kinetic_energy(water_basis);
    EXPECT_TRUE(water_expected.isApprox(water_actual, 1e-8));

    Eigen::MatrixXd ethene_expected(ethene_basis.size(), ethene_basis.size());
    ethene_expected <<
               15.891122,  -0.085889994,             0,             0,             0, -8.5290384e-05, -0.0083193283,  0.0070135536, -0.0070604917, 0.00033749263,  -0.010371908,  -0.010377125, -0.0025727365, -0.0025727458,
            -0.085889994,    0.47224999,             0,             0,             0, -0.0083193283,    0.05819253,   -0.12076679,    0.12157502, -0.0058113055,    0.10825678,    0.10821026,   -0.01419343,  -0.014193388,
                       0,             0,     1.4777281,             0,             0, -0.0070135536,    0.12076679,   -0.12511262,    0.21842205,  -0.010440609,   0.061229979,   -0.25564952,  0.0038103234,   0.001321334,
                       0,             0,             0,     1.4777281,             0,  0.0070604917,   -0.12157502,    0.21842205,   -0.12802647,   0.010510482,    0.25547849,  -0.063062256, -0.0013537952, -0.0038123201,
                       0,             0,             0,             0,     1.4777281, -0.00033749263,  0.0058113055,  -0.010440609,   0.010510482,   0.091354958,  -0.023943286,   0.014558994, -2.4709453e-05, 0.00027166036,
          -8.5290384e-05, -0.0083193283, -0.0070135536,  0.0070604917, -0.00033749263,     15.891122,  -0.085889994,             0,             0,             0, -0.0024128851,  -0.002438028, -0.0077410349, -0.0077407962,
           -0.0083193283,    0.05819253,    0.12076679,   -0.12157502,  0.0058113055,  -0.085889994,    0.47224999,             0,             0,             0,  -0.014762137,  -0.014692407,    0.12675347,    0.12675486,
            0.0070135536,   -0.12076679,   -0.12511262,    0.21842205,  -0.010440609,             0,             0,     1.4777281,             0,             0, -0.0002949934, -0.0012971915,    0.28780429,  -0.077987093,
           -0.0070604917,    0.12157502,    0.21842205,   -0.12802647,   0.010510482,             0,             0,             0,     1.4777281,             0, 0.00085417345, 0.00044826395,   0.075027423,   -0.28625192,
           0.00033749263, -0.0058113055,  -0.010440609,   0.010510482,   0.091354958,             0,             0,             0,             0,     1.4777281, -6.1442726e-05, 9.4719675e-06,  -0.016727605,   0.026825112,
            -0.010371908,    0.10825678,   0.061229979,    0.25547849,  -0.023943286, -0.0024128851,  -0.014762137, -0.0002949934, 0.00085417345, -6.1442726e-05,    0.76003188, -0.0076419431,  -0.011052475, -0.0063714232,
            -0.010377125,    0.10821026,   -0.25564952,  -0.063062256,   0.014558994,  -0.002438028,  -0.014692407, -0.0012971915, 0.00044826395, 9.4719675e-06, -0.0076419431,    0.76003188, -0.0064291459,  -0.011123053,
           -0.0025727365,   -0.01419343,  0.0038103234, -0.0013537952, -2.4709453e-05, -0.0077410349,    0.12675347,    0.28780429,   0.075027423,  -0.016727605,  -0.011052475, -0.0064291459,    0.76003188, -0.0058082006,
           -0.0025727458,  -0.014193388,   0.001321334, -0.0038123201, 0.00027166036, -0.0077407962,    0.12675486,  -0.077987093,   -0.28625192,   0.026825112, -0.0063714232,  -0.011123053, -0.0058082006,    0.76003188;

    auto ethene_actual = kinetic_energy(ethene_basis);
    EXPECT_TRUE(ethene_expected.isApprox(ethene_actual, 1e-7));
}

TEST(OneElectronIntegralTest, potential_energy_matrix_test) {
    Eigen::MatrixXd water_expected(water_basis.size(), water_basis.size());
    water_expected <<    -61.724045,    -7.4447756,             0,  -0.019000303,             0,    -1.7467216,    -1.7467216,
                         -7.4447756,    -10.142872,             0,   -0.22349655,             0,     -3.869445,     -3.869445,
                                  0,             0,    -10.142626,             0,             0,    -2.2547736,     2.2547736,
                       -0.019000303,   -0.22349655,             0,    -10.079964,             0,    -1.8181924,    -1.8181924,
                                  0,             0,             0,             0,    -9.9863261,             0,             0,
                         -1.7467216,     -3.869445,    -2.2547736,    -1.8181924,             0,    -5.8364114,    -1.6165044,
                         -1.7467216,     -3.869445,     2.2547736,    -1.8181924,             0,    -1.6165044,    -5.8364114;

    auto water_actual = potential_energy(water_basis, water);
    EXPECT_TRUE(water_expected.isApprox(water_actual, 1e-8));

    Eigen::MatrixXd ethene_expected(ethene_basis.size(), ethene_basis.size());
    ethene_expected <<
             -37.51299,    -5.1509695,  -0.040882385,   0.041371281, -0.0019655294, -2.6748192e-05,   -0.91889667,     1.0750219,    -1.0822058,   0.051729942,    -1.3222636,     -1.321933,   -0.14500132,   -0.14500458,
            -5.1509695,    -9.0014883,   -0.46893857,    0.47428687,  -0.022548993,   -0.92034389,    -3.6104171,     2.7302431,    -2.7479323,    0.13135471,    -3.7371281,    -3.7403085,   -0.94830396,   -0.94850738,
          -0.040882385,   -0.46893857,    -9.0881502,    0.20298809, -0.0025903928,    -1.0767026,    -2.7337545,    0.70763102,    -2.6123635,    0.12602784,   -0.89700399,     3.0921709,    -1.0877242,   -0.46596317,
           0.041371281,    0.47428687,    0.20298809,    -9.0917949,   0.026349197,     1.0839119,     2.7528029,    -2.6133512,    0.74302106,   -0.12308521,    -3.0862372,    0.92358643,    0.47513873,     1.0894442,
         -0.0019655294,  -0.022548993, -0.0025903928,   0.026349197,    -8.6934946,  -0.051811158,   -0.13157725,    0.12606757,   -0.12307769,    -1.8502005,    0.29514633,    -0.1895283, -0.00036927079,  -0.074421378,
        -2.6748192e-05,   -0.92034389,    -1.0767026,     1.0839119,  -0.051811158,      -37.5497,    -5.1600869,   0.039858617,  -0.040055103,   0.001914708,   -0.12424717,   -0.12745422,    -1.5049495,    -1.5049736,
           -0.91889667,    -3.6104171,    -2.7337545,     2.7528029,   -0.13157725,    -5.1600869,    -9.0288717,     0.4621907,   -0.46446719,   0.022202589,    -0.8407587,   -0.85566234,    -4.0499949,    -4.0504289,
             1.0750219,     2.7302431,    0.70763102,    -2.6133512,    0.12606757,   0.039858617,     0.4621907,    -9.1205463,    0.19055667, -0.0011641031,    0.41582996,    0.99355562,    -3.2573362,     1.0601077,
            -1.0822058,    -2.7479323,    -2.6123635,    0.74302106,   -0.12307769,  -0.040055103,   -0.46446719,    0.19055667,    -9.1213505,   0.027409955,   -0.97715736,   -0.42334375,    -1.0256394,     3.2385655,
           0.051729942,    0.13135471,    0.12602784,   -0.12308521,    -1.8502005,   0.001914708,   0.022202589, -0.0011641031,   0.027409955,    -8.7143431,   0.067371159, -0.00054326619,     0.2041281,   -0.30993061,
            -1.3222636,    -3.7371281,   -0.89700399,    -3.0862372,    0.29514633,   -0.12424717,    -0.8407587,    0.41582996,   -0.97715736,   0.067371159,    -6.2605175,     -1.034323,   -0.37804898,   -0.13929903,
             -1.321933,    -3.7403085,     3.0921709,    0.92358643,    -0.1895283,   -0.12745422,   -0.85566234,    0.99355562,   -0.42334375, -0.00054326619,     -1.034323,    -6.2691746,   -0.14134712,   -0.38886621,
           -0.14500132,   -0.94830396,    -1.0877242,    0.47513873, -0.00036927079,    -1.5049495,    -4.0499949,    -3.2573362,    -1.0256394,     0.2041281,   -0.37804898,   -0.14134712,    -6.4497798,     -1.186581,
           -0.14500458,   -0.94850738,   -0.46596317,     1.0894442,  -0.074421378,    -1.5049736,    -4.0504289,     1.0601077,     3.2385655,   -0.30993061,   -0.13929903,   -0.38886621,     -1.186581,    -6.4508393;

    auto ethene_actual = potential_energy(ethene_basis, ethene);
    EXPECT_TRUE(ethene_expected.isApprox(ethene_actual, 1e-8));
}

class OneElectronIntegralDeathTest: public IntegralTestData{};

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
