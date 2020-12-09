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

#include <unordered_map>

#include "CommonTestingData.hpp"
#include "IntegralTestData.hpp"
#include "TwoElectronTestData.hpp"

#include "TwoElectronIntegralDetail.hpp"

struct TwoElectronIntegralIndexTuple {
    uint16_t mu, nu, lambda, sigma;

    bool operator==(const TwoElectronIntegralIndexTuple& other) const = default;
};

struct IndexTupleHash {
    auto operator()(const TwoElectronIntegralIndexTuple& index_tuple) const{
        auto [mu, nu, lambda, sigma] = index_tuple;
        std::uint64_t hash_input =
                  (static_cast<uint64_t>(mu))
                | (static_cast<uint64_t>(nu))     << 16u
                | (static_cast<uint64_t>(lambda)) << 32u
                | (static_cast<uint64_t>(sigma))  << 48u;
        return std::hash<uint64_t>{}(hash_input);
    }
};

using TwoElectronIntegralVector = std::vector<TwoElectronIntegral>;
using TwoElectronIntegralMap = std::unordered_map<TwoElectronIntegralIndexTuple, double, IndexTupleHash>;

TwoElectronIntegralMap integral_map(const TwoElectronIntegralVector& integral_vec) {
    TwoElectronIntegralMap result;
    for (auto [mu, nu, lambda, sigma, val]: integral_vec) {
        result[TwoElectronIntegralIndexTuple{mu, nu, lambda, sigma}] = val;
    }
    return result;
}

double integral_or_zero(const TwoElectronIntegralMap& integral_map, uint16_t mu, uint16_t nu, uint16_t lambda, uint16_t sigma) {
    TwoElectronIntegralIndexTuple index_tuple{mu, nu, lambda, sigma};
    if (integral_map.contains(index_tuple))
        return integral_map.at(index_tuple);
    else
        return 0;
}

void test_two_electron_integrals(const TwoElectronIntegralVector& test, const TwoElectronIntegralVector& reference, uint16_t n) {
    ASSERT_EQ((n*n*n*n + 2*n*n*n + 3*n*n + 2*n)/8, reference.size());
    auto test_map = integral_map(test);
    for (auto [mu, nu, lambda, sigma, val]: reference) {
        EXPECT_NEAR(integral_or_zero(test_map, mu, nu, lambda, sigma), val, 1e-10);
    }
}

TEST(TwoElectronIntegralTest, water_integrals) {
    test_two_electron_integrals(calculate_two_electron_integrals(water_basis), water_ref, water_basis.size());
}

TEST(TwoElectronIntegralTest, ethene_integrals) {
    test_two_electron_integrals(calculate_two_electron_integrals(ethene_basis), ethene_ref, ethene_basis.size());
}

class TwoElectronIntegralDeathTest: public IntegralTestData{};

TEST_F(TwoElectronIntegralDeathTest, integral_ssss) {
    EXPECT_DEBUG_DEATH(coulomb_ssss(p, s, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(s, p, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(p, p, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(s, s, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(p, s, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(s, p, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(p, p, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(s, s, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(p, s, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(s, p, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(p, p, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(s, s, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(p, s, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(s, p, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ssss(p, p, p, p), "");
}

TEST_F(TwoElectronIntegralDeathTest, integral_psss) {
    EXPECT_DEBUG_DEATH(coulomb_psss(s, s, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(s, p, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(p, p, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(s, s, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(p, s, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(s, p, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(p, p, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(s, s, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(p, s, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(s, p, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(p, p, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(s, s, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(p, s, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(s, p, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psss(p, p, p, p), "");
}

TEST_F(TwoElectronIntegralDeathTest, integral_ppss) {
    EXPECT_DEBUG_DEATH(coulomb_ppss(s, s, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(p, s, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(s, p, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(s, s, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(p, s, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(s, p, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(p, p, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(s, s, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(p, s, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(s, p, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(p, p, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(s, s, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(p, s, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(s, p, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppss(p, p, p, p), "");
}

TEST_F(TwoElectronIntegralDeathTest, integral_psps) {
    EXPECT_DEBUG_DEATH(coulomb_psps(s, s, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(p, s, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(s, p, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(p, p, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(s, s, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(s, p, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(p, p, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(s, s, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(p, s, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(s, p, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(p, p, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(s, s, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(p, s, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(s, p, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_psps(p, p, p, p), "");
}

TEST_F(TwoElectronIntegralDeathTest, integral_ppps) {
    EXPECT_DEBUG_DEATH(coulomb_ppps(s, s, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(p, s, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(s, p, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(p, p, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(s, s, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(p, s, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(s, p, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(s, s, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(p, s, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(s, p, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(p, p, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(s, s, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(p, s, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(s, p, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_ppps(p, p, p, p), "");
}

TEST_F(TwoElectronIntegralDeathTest, integral_pppp) {
    EXPECT_DEBUG_DEATH(coulomb_pppp(s, s, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(p, s, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(s, p, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(p, p, s, s), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(s, s, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(p, s, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(s, p, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(p, p, p, s), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(s, s, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(p, s, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(s, p, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(p, p, s, p), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(s, s, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(p, s, p, p), "");
    EXPECT_DEBUG_DEATH(coulomb_pppp(s, p, p, p), "");
}
