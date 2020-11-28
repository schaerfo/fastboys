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

#ifndef COMMONTESTINGDATA_HPP
#define COMMONTESTINGDATA_HPP

#include <Eigen/Dense>
#include "Molecule.hpp"

const Eigen::Vector3d water_o_pos{0, 0, 0};
const Eigen::Vector3d water_h1_pos{1.430541, 1.1072565, 0};
const Eigen::Vector3d water_h2_pos{-1.430541, 1.1072565, 0};

const Molecule water({Atom{"O", water_o_pos},
                      Atom{"H", water_h1_pos},
                      Atom{"H", water_h2_pos}});

using PrimFunc = BasisFunction::PrimitiveFunction;

const BasisSet water_basis{
        BasisFunction{water_o_pos, {PrimFunc{130.7093214, 4.251943277890349}, PrimFunc{23.80886605, 4.112294424303829}, PrimFunc{6.443608313, 1.2816225514654456}}, BasisFunction::Type::s},
        BasisFunction{water_o_pos, {PrimFunc{5.033151319, -0.2394130049495816}, PrimFunc{1.169596125, 0.3202342354447327}, PrimFunc{0.38038896, 0.24168555456024995}}, BasisFunction::Type::s},
        BasisFunction{water_o_pos, {PrimFunc{5.033151319, 1.6754501961028012}, PrimFunc{1.169596125, 1.0535680440009998}, PrimFunc{0.38038896, 0.1669028790864184}}, BasisFunction::Type::px},
        BasisFunction{water_o_pos, {PrimFunc{5.033151319, 1.6754501961028012}, PrimFunc{1.169596125, 1.0535680440009998}, PrimFunc{0.38038896, 0.1669028790864184}}, BasisFunction::Type::py},
        BasisFunction{water_o_pos, {PrimFunc{5.033151319, 1.6754501961028012}, PrimFunc{1.169596125, 1.0535680440009998}, PrimFunc{0.38038896, 0.1669028790864184}}, BasisFunction::Type::pz},
        BasisFunction{water_h1_pos, {PrimFunc{3.425250914, 0.2769343550887434}, PrimFunc{0.6239137298, 0.26783885161885207}, PrimFunc{0.168855404, 0.08347367113276243}}, BasisFunction::Type::s},
        BasisFunction{water_h2_pos, {PrimFunc{3.425250914, 0.2769343550887434}, PrimFunc{0.6239137298, 0.26783885161885207}, PrimFunc{0.168855404, 0.08347367113276243}}, BasisFunction::Type::s}
};

const Eigen::Vector3d ethene_c1_pos{-6.0221637, 3.4502706, 0.0054999};
const Eigen::Vector3d ethene_c2_pos{-4.2401772, 1.6563582, 0.0912492};
const Eigen::Vector3d ethene_h1_pos{-5.5459404, 5.4372843, -0.1807218};
const Eigen::Vector3d ethene_h2_pos{-8.0114265, 2.9595699, 0.1187865};
const Eigen::Vector3d ethene_h3_pos{-2.3557527, 2.147607, -0.0182763};
const Eigen::Vector3d ethene_h4_pos{-4.7507985, -0.2178792, 0.2668869};

const Molecule ethene({Atom("C", ethene_c1_pos),
                       Atom("C", ethene_c2_pos),
                       Atom("H", ethene_h1_pos),
                       Atom("H", ethene_h2_pos),
                       Atom("H", ethene_h3_pos),
                       Atom("H", ethene_h4_pos)});

const BasisSet ethene_basis{
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

#endif //COMMONTESTINGDATA_HPP
