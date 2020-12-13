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

#include <thread>

#include "CommonTestingData.hpp"
#include "TwoElectronTestData.hpp"
#include "FockMatrix.hpp"

TEST(FockMatrixTest, water) {
    Eigen::MatrixXd p(water_basis.size(), water_basis.size());
    auto f_expected = p;
    p  <<    2.09606032244,     -0.80675351078,  6.42074641652e-17,     -0.14627973355,                  0,     0.069881056854,     0.069881056854,
            -0.80675351078,       3.1787064872, -1.68597705613e-16,     0.532519573843,  3.01269625005e-32,     -1.75952105849,     -1.75952105849,
         6.42074641652e-17, -1.68597705613e-16,      1.74060273453,  1.26314619506e-17, -8.00913031675e-16,     -1.63214918838,      1.63214918838,
            -0.14627973355,     0.532519573843,  1.26314619506e-17,       2.2670739999, -2.53543020578e-31,    -0.856011180761,    -0.856011180761,
                         0,  3.01269625005e-32, -8.00913031675e-16, -2.53543020578e-31,      2.00000000004, -2.75114524036e-16,  2.75114524036e-16,
            0.069881056854,     -1.75952105849,     -1.63214918838,    -0.856011180761, -2.75114524036e-16,      2.72521351247,    -0.335692827957,
            0.069881056854,     -1.75952105849,      1.63214918838,    -0.856011180761,  2.75114524036e-16,    -0.335692827957,      2.72521351247;

    f_expected <<      12.394404814,       2.6240454262,  3.46944695195e-18,    0.0969207662784,                  0,     0.683467028183,     0.683467028183,
                       2.6240454262,      7.13089123599,  5.55111512313e-17,     0.241016570089,  3.08148791102e-33,      3.03127860372,      3.03127860372,
                  3.46944695195e-18,  5.55111512313e-17,      7.46098746283, -2.77555756156e-17,  3.02287609033e-16,      1.90595984446,     -1.90595984446,
                    0.0969207662784,     0.241016570089, -2.77555756156e-17,      7.01380694579,  8.55112895308e-32,      1.46536929024,      1.46536929024,
                                  0,  3.08148791102e-33,  3.02287609033e-16,  8.55112895308e-32,      6.89904894007,  1.10467866121e-16, -1.10467866121e-16,
                     0.683467028183,      3.03127860372,      1.90595984446,      1.46536929024,  1.10467866121e-16,      4.73429237209,      1.34017800307,
                     0.683467028183,      3.03127860372,     -1.90595984446,      1.46536929024, -1.10467866121e-16,      1.34017800307,      4.73429237209;

    auto f_actual = electron_repulsion_matrix(water_ref, p, std::thread::hardware_concurrency());
    ASSERT_EQ(f_actual.cols(), f_expected.cols());
    ASSERT_EQ(f_actual.rows(), f_expected.rows());
    EXPECT_TRUE(f_expected.isApprox(f_actual, 1e-11));
}

TEST(FockMatrixTest, ethene) {
    Eigen::MatrixXd p(ethene_basis.size(), ethene_basis.size());
    auto f_expected = p;
    p  <<      2.1283087963,     -1.09752208423,   -0.0520894781394,    0.0525284661922,  -0.00263058003476,   -0.0236829427155,     0.255879471539,   -0.0440877794089,    0.0445388861344,  -0.00242000720436,     0.152368632692,     0.154549868044,    0.0279581250425,    0.0276201321846,
             -1.09752208423,      5.42328560712,     0.454241187779,    -0.473086316461,    0.0219819882056,     0.257426519826,      -2.3409498713,      1.60796987927,     -1.62081606317,    0.0774100792186,       -2.138658849,     -2.14276913822,    0.0453761868128,    0.0468000161338,
           -0.0520894781394,     0.454241187779,      2.90422796022,    -0.235902789178,   -0.0382577000243,    0.0519986505074,     -1.60925835703,     0.678193858827,     -1.49808200456,    0.0532743164867,    -0.589969475283,      2.28218148076,     0.237485502237,    -0.204371826853,
            0.0525284661922,    -0.473086316461,    -0.235902789178,      2.92237015492,    -0.127590233562,   -0.0552304833514,      1.62474059959,     -1.50107184539,     0.692537805064,    -0.114033349911,     -2.28276095155,     0.606337098724,     0.199935830188,    -0.236798577002,
          -0.00263058003476,    0.0219819882056,   -0.0382577000243,    -0.127590233562,       1.3241222721,   0.00278107458656,   -0.0778308460523,    0.0537122644427,     -0.11366990901,     -1.30430869585,     0.215209387338,    -0.133670829037,   -0.0254789037197,    0.0272818145558,
           -0.0236829427155,     0.257426519826,    0.0519986505074,   -0.0552304833514,   0.00278107458656,      2.13518762989,     -1.12803940156,    0.0431141719935,   -0.0421730419173,   0.00203391061622,    0.0260327122602,    0.0270869810508,     0.175620493985,     0.174209929329,
             0.255879471539,      -2.3409498713,     -1.60925835703,      1.62474059959,   -0.0778308460523,     -1.12803940156,      5.90586318165,    -0.359933441222,     0.361667988131,   -0.0172557303585,    0.0339036426594,    0.0459906477045,     -2.45260769327,     -2.45376113953,
           -0.0440877794089,      1.60796987927,     0.678193858827,     -1.50107184539,    0.0537122644427,    0.0431141719935,    -0.359933441222,      3.07064153676,    -0.107466987753,   -0.0547561978091,     0.210456589953,     -0.22888106222,     -2.51105395453,     0.736013168926,
            0.0445388861344,     -1.62081606317,     -1.49808200456,     0.692537805064,     -0.11366990901,   -0.0421730419173,     0.361667988131,    -0.107466987753,      3.06237863453,     -0.14274077565,     0.223458459057,    -0.201855086517,    -0.709137622028,      2.49612344615,
          -0.00242000720436,    0.0774100792186,    0.0532743164867,    -0.114033349911,     -1.30430869585,   0.00203391061622,   -0.0172557303585,   -0.0547561978091,     -0.14274077565,      1.32587407018,   -0.0265792765749,    0.0252995715772,     0.150613609175,    -0.236095672156,
             0.152368632692,       -2.138658849,    -0.589969475283,     -2.28276095155,     0.215209387338,    0.0260327122602,    0.0339036426594,     0.210456589953,     0.223458459057,   -0.0265792765749,      3.17909211287,   -0.0025006011769,   -0.0567410161332,   -0.0812431949623,
             0.154549868044,     -2.14276913822,      2.28218148076,     0.606337098724,    -0.133670829037,    0.0270869810508,    0.0459906477045,     -0.22888106222,    -0.201855086517,    0.0252995715772,   -0.0025006011769,      3.17743850076,   -0.0838772975326,   -0.0590030999718,
            0.0279581250425,    0.0453761868128,     0.237485502237,     0.199935830188,   -0.0254789037197,     0.175620493985,     -2.45260769327,     -2.51105395453,    -0.709137622028,     0.150613609175,   -0.0567410161332,   -0.0838772975326,      3.57427183987,   -0.0335192555466,
            0.0276201321846,    0.0468000161338,    -0.204371826853,    -0.236798577002,    0.0272818145558,     0.174209929329,     -2.45376113953,     0.736013168926,      2.49612344615,    -0.236095672156,   -0.0812431949623,   -0.0590030999718,   -0.0335192555466,      3.57491423777;

    f_expected <<     11.038969771,      2.55473624292,      0.04869990978,   -0.0484636907707,   0.00237726749894,    0.0100739982866,     0.485303143426,    -0.568472931631,     0.572359237946,   -0.0272910918033,     0.697794463494,     0.697171033831,    0.0730133682798,    0.0730086078824,
                     2.55473624292,      7.44773075983,      0.40929546623,     -0.41346796963,     0.019675823996,      0.48581146123,      3.08780522623,     -2.33869374921,      2.35391230076,    -0.112485287752,      3.25117168683,      3.25372137204,     0.831449405592,     0.831605043876,
                     0.04869990978,      0.40929546623,      7.58565504489,    -0.155700201161,  0.000736348193165,     0.567817040882,      2.34095843453,    -0.602057611613,      2.30542399886,    -0.106416972278,     0.804122961967,     -2.76878777925,     0.964702017017,     0.388044420345,
                  -0.0484636907707,     -0.41346796963,    -0.155700201161,      7.58768240939,   -0.0230429214334,     -0.57097747975,     -2.35675208155,      2.30575802083,    -0.632088221272,     0.119658534666,      2.76307444088,    -0.827276977764,     -0.39617707593,    -0.965982794693,
                  0.00237726749894,     0.019675823996,  0.000736348193165,   -0.0230429214334,      7.24878762721,    0.0272616028272,      0.11263586637,    -0.106431309173,     0.119570930228,      1.78581828306,    -0.264224936228,     0.169745258738,  -0.00177819124352,    0.0668955798038,
                   0.0100739982866,      0.48581146123,     0.567817040882,     -0.57097747975,    0.0272616028272,      11.0804198698,      2.56313970554,   -0.0445695426681,     0.044268238309,  -0.00212458903527,    0.0622644688319,    0.0637931539281,     0.786508673886,     0.786955637193,
                    0.485303143426,      3.08780522623,      2.34095843453,     -2.35675208155,      0.11263586637,      2.56313970554,      7.47124864479,    -0.404094350728,     0.405947876271,   -0.0193975904072,     0.738863698958,      0.75226074372,      3.50375602136,      3.50429643341,
                   -0.568472931631,     -2.33869374921,    -0.602057611613,      2.30575802083,    -0.106431309173,   -0.0445695426681,    -0.404094350728,      7.60944052287,    -0.147242643588, -0.000221143229694,    -0.347182217786,    -0.883341743914,      2.89560611559,     -0.94305034356,
                    0.572359237946,      2.35391230076,      2.30542399886,    -0.632088221272,     0.119570930228,     0.044268238309,     0.405947876271,    -0.147242643588,      7.60991308219,   -0.0237760923317,     0.868144239908,     0.353950090361,     0.912297973861,     -2.87894865871,
                  -0.0272910918033,    -0.112485287752,    -0.106416972278,     0.119658534666,      1.78581828306,  -0.00212458903527,   -0.0193975904072, -0.000221143229694,   -0.0237760923317,      7.26497630347,   -0.0607027333705,     0.002356814397,    -0.181518143013,     0.275567064278,
                    0.697794463494,      3.25117168683,     0.804122961967,      2.76307444088,    -0.264224936228,    0.0622644688319,     0.738863698958,    -0.347182217786,     0.868144239908,   -0.0607027333705,      5.37846173809,     0.916557960307,     0.341786868436,     0.126840172594,
                    0.697171033831,      3.25372137204,     -2.76878777925,    -0.827276977764,     0.169745258738,    0.0637931539281,      0.75226074372,    -0.883341743914,     0.353950090361,     0.002356814397,     0.916557960307,      5.38617664512,     0.128661659976,     0.351611830387,
                   0.0730133682798,     0.831449405592,     0.964702017017,     -0.39617707593,  -0.00177819124352,     0.786508673886,      3.50375602136,      2.89560611559,     0.912297973861,    -0.181518143013,     0.341786868436,     0.128661659976,      5.54084791896,      1.04297404785,
                   0.0730086078824,     0.831605043876,     0.388044420345,    -0.965982794693,    0.0668955798038,     0.786955637193,      3.50429643341,     -0.94305034356,     -2.87894865871,     0.275567064278,     0.126840172594,     0.351611830387,      1.04297404785,      5.54203237129;

    auto f_actual = electron_repulsion_matrix(ethene_ref, p, std::thread::hardware_concurrency());
    ASSERT_EQ(f_actual.cols(), f_expected.cols());
    ASSERT_EQ(f_actual.rows(), f_expected.rows());
    EXPECT_TRUE(f_expected.isApprox(f_actual, 1e-11));
}
