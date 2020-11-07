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

#include "Helpers.hpp"

class HelpersTest: public ::testing::Test {
protected:
    void SetUp() override {
        s = px = py = pz = BasisFunction{{0, 0, 0}, {}, BasisFunction::Type::s};
        px.type = BasisFunction::Type::px;
        py.type = BasisFunction::Type::py;
        pz.type = BasisFunction::Type::pz;
    }

    BasisFunction s, px, py, pz;
};

TEST_F(HelpersTest, is_s_orbital_test) {
    EXPECT_TRUE(helpers::is_s_orbital(s.type));
    EXPECT_FALSE(helpers::is_s_orbital(px.type));
    EXPECT_FALSE(helpers::is_s_orbital(py.type));
    EXPECT_FALSE(helpers::is_s_orbital(pz.type));
}

TEST_F(HelpersTest, is_p_orbital_test) {
    EXPECT_FALSE(helpers::is_p_orbital(s.type));
    EXPECT_TRUE(helpers::is_p_orbital(px.type));
    EXPECT_TRUE(helpers::is_p_orbital(py.type));
    EXPECT_TRUE(helpers::is_p_orbital(pz.type));
}

TEST_F(HelpersTest, vector_component_test) {
    px.center = {0, 1, 2};
    EXPECT_EQ(helpers::vector_component(px.center, BasisFunction::Type::px), 0);
    EXPECT_EQ(helpers::vector_component(px.center, BasisFunction::Type::py), 1);
    EXPECT_EQ(helpers::vector_component(px.center, BasisFunction::Type::pz), 2);
}

TEST_F(HelpersTest, component_difference_test) {
    px.center = {1, 2, 3};
    EXPECT_EQ(helpers::component_difference(px.center, py.center, BasisFunction::Type::px), 1);
    EXPECT_EQ(helpers::component_difference(px.center, py.center, BasisFunction::Type::py), 2);
    EXPECT_EQ(helpers::component_difference(px.center, py.center, BasisFunction::Type::pz), 3);

    EXPECT_EQ(helpers::component_difference(py.center, px.center, BasisFunction::Type::px), -1);
    EXPECT_EQ(helpers::component_difference(py.center, px.center, BasisFunction::Type::py), -2);
    EXPECT_EQ(helpers::component_difference(py.center, px.center, BasisFunction::Type::pz), -3);
}
