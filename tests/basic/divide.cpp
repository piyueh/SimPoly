/**
 * \file divide.cpp
 * \brief Unit tests for polynomial division.
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
