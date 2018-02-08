/**
 * \file operators2.cpp
 * \brief Unit tests for polynomial basic operations.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-29
 */


# include <random>

# include <gtest/gtest.h>

# include "basic.h"
# include "exceptions.h"


extern std::default_random_engine generator;

using namespace simpoly;

# ifndef NDEBUG
TEST(PolynomialMult, ZeroLengthTests)
{
    std::valarray<double> p1(0), p2(1);
    
    try { basic::multiply(p1, p2); }
    catch (simpoly::exceptions::ZeroCoeffsLength e) {};
    
    try { basic::multiply(p2, p1); }
    catch (simpoly::exceptions::ZeroCoeffsLength e) {};
}
# endif

TEST(PolynomialMult, FixedPolynomials1)
{
    std::valarray<double> p1({0.0, 0.0, 3.0});
    std::valarray<double> p2({7.0, -5.0, 4.0});
    std::valarray<double> expect({0.0, 0.0, 21.0, -15.0, 12.0});
    std::valarray<double> actual;
    
    actual = basic::multiply(p1, p2);
    ASSERT_EQ(5, actual.size());
    for(unsigned i=0; i<5; ++i) ASSERT_NEAR(expect[i], actual[i], 1e-12);
    
    actual = basic::multiply(p2, p1);
    ASSERT_EQ(5, actual.size());
    for(unsigned i=0; i<5; ++i) ASSERT_NEAR(expect[i], actual[i], 1e-12);
}

TEST(PolynomialMult, FixedPolynomials2)
{
    std::valarray<double> p1({-5.0, 4.0});
    std::valarray<double> p2({-6.0, 3.0, 2.0});
    std::valarray<double> expect({30.0, -39.0, 2.0, 8.0});
    std::valarray<double> actual;
    
    actual = basic::multiply(p1, p2);
    ASSERT_EQ(4, actual.size());
    for(unsigned i=0; i<4; ++i) ASSERT_NEAR(expect[i], actual[i], 1e-12);
    
    actual = basic::multiply(p2, p1);
    ASSERT_EQ(4, actual.size());
    for(unsigned i=0; i<4; ++i) ASSERT_NEAR(expect[i], actual[i], 1e-12);
}

TEST(PolynomialMult, FixedPolynomials3)
{
    std::valarray<double> p1({5.0, 0.0, 1.0});
    std::valarray<double> p2({9.0, -19.0, 1.0});
    std::valarray<double> expect({45.0, -95.0, 14.0, -19.0, 1.0});
    std::valarray<double> actual;
    
    actual = basic::multiply(p1, p2);
    ASSERT_EQ(5, actual.size());
    for(unsigned i=0; i<5; ++i) ASSERT_NEAR(expect[i], actual[i], 1e-12);
    
    actual = basic::multiply(p2, p1);
    ASSERT_EQ(5, actual.size());
    for(unsigned i=0; i<5; ++i) ASSERT_NEAR(expect[i], actual[i], 1e-12);
}

TEST(PolynomialMult, RandomPolynomials)
{
    std::uniform_int_distribution<int> dist1(1, 20);
    std::uniform_real_distribution<double> dist2(-5.0, 5.0);
    
    int degree1 = dist1(generator),
        degree2 = dist1(generator),
        result_degree = degree1 + degree2;
    
    std::valarray<double> p1(degree1+1);
    std::valarray<double> p2(degree2+1);
    std::valarray<double> actual;
    
    for(unsigned i=0; i<degree1+1; ++i) p1[i] = dist2(generator);
    for(unsigned i=0; i<degree2+1; ++i) p2[i] = dist2(generator);
    
    actual = basic::multiply(p1, p2);
    ASSERT_EQ(result_degree+1, actual.size());
    for(unsigned i=0; i<result_degree+1; ++i)
    {
        double value = 0;
        for(unsigned j=0; j<degree1+1; ++j)
        {
            if (i < j) break;
            if ((i - j) > degree2) continue;
            value += (p1[j] * p2[i-j]);
        }
        ASSERT_NEAR(value, actual[i], 1e-12);
    }
    
    actual = basic::multiply(p2, p1);
    ASSERT_EQ(result_degree+1, actual.size());
    for(unsigned i=0; i<result_degree+1; ++i)
    {
        double value = 0;
        for(unsigned j=0; j<degree1+1; ++j)
        {
            if (i < j) break;
            if ((i - j) > degree2) continue;
            value += (p1[j] * p2[i-j]);
        }
        ASSERT_NEAR(value, actual[i], 1e-12);
    }
}


# ifndef NDEBUG
TEST(PolynomialDivide, ZeroLengthTests)
{
    std::valarray<double> p1(0), p2(1);
    
    try { basic::divide(p1, p2); }
    catch (exceptions::ZeroCoeffsLength e) {};
    
    try { basic::divide(p2, p1); }
    catch (exceptions::ZeroCoeffsLength e) {};
}
# endif

TEST(PolynomialDivide, LongerSecondPoly)
{
    std::valarray<double> p1({1.0, 2.0, 3.0}), p2({1.0, 2.0, 3.0, 4.0});
    
    // test the version without remainder
    std::valarray<double> q = basic::divide(p1, p2);
    ASSERT_EQ(1, q.size());
    ASSERT_DOUBLE_EQ(0.0, q[0]);
    
    // test the version with remainder
    std::valarray<double> r;
    q = basic::divide(p1, p2, r);
    ASSERT_EQ(1, q.size());
    ASSERT_DOUBLE_EQ(0.0, q[0]);
    for(unsigned i=0; i<3; ++i) ASSERT_DOUBLE_EQ(p1[i], r[i]);
}

TEST(PolynomialDivide, FixPolynomial1)
{
    std::valarray<double> p1({-10.0, -9.0, 1.0}), p2({1.0, 1.0});
    std::valarray<double> q, r;
    
    q = basic::divide(p1, p2, r);
    
    ASSERT_EQ(2, q.size());
    ASSERT_DOUBLE_EQ(-10.0, q[0]);
    ASSERT_DOUBLE_EQ(1.0, q[1]);
    
    ASSERT_EQ(1, r.size());
    ASSERT_DOUBLE_EQ(0.0, r[0]);
}

TEST(PolynomialDivide, FixPolynomial2)
{
    std::valarray<double> p1({-3.0, 10.0, -5.0, 3.0}), p2({1.0, 3.0});
    std::valarray<double> q, r;
    
    q = basic::divide(p1, p2, r);
    
    ASSERT_EQ(3, q.size());
    ASSERT_DOUBLE_EQ(4.0, q[0]);
    ASSERT_DOUBLE_EQ(-2.0, q[1]);
    ASSERT_DOUBLE_EQ(1.0, q[2]);
    
    ASSERT_EQ(1, r.size());
    ASSERT_DOUBLE_EQ(-7.0, r[0]);
}

TEST(PolynomialDivide, FixPolynomial3)
{
    std::valarray<double> p1({1.0, 2.0, 0.0, 3.0, 4.0}), p2({2.0, 1.0, 1.0});
    std::valarray<double> q, r;
    
    q = basic::divide(p1, p2, r);
    
    ASSERT_EQ(3, q.size());
    ASSERT_DOUBLE_EQ(-7.0, q[0]);
    ASSERT_DOUBLE_EQ(-1.0, q[1]);
    ASSERT_DOUBLE_EQ(4.0, q[2]);
    
    ASSERT_EQ(2, r.size());
    ASSERT_DOUBLE_EQ(15.0, r[0]);
    ASSERT_DOUBLE_EQ(11.0, r[1]);
}

TEST(PolynomialDivide, RandomPolynomial)
{
    std::uniform_int_distribution<int> dist1(1, 20);
    int len1 = dist1(generator);
    
    std::uniform_int_distribution<int> dist2(1, len1);
    int len2 = dist2(generator);
    
    std::valarray<double> p1(len1), p2(len2);
    std::uniform_real_distribution<double> rdist(1, 2);
    std::uniform_int_distribution<int> bdist(0, 1);
    for(unsigned i=0; i<len1; ++i) 
        p1[i] = rdist(generator)*((bdist(generator))? 1.0 : -1.0);
    for(unsigned i=0; i<len2; ++i) 
        p2[i] = rdist(generator)*((bdist(generator))? 1.0 : -1.0);
    
    std::valarray<double> q, r, expect;
    
    q = basic::divide(p1, p2, r);
    
    expect = basic::add(basic::multiply(p2, q), r);
    ASSERT_EQ(len1, expect.size());
    
    for(unsigned i=0; i<len1; ++i) ASSERT_NEAR(p1[i], expect[i], 1e-8);
}
