/**
 * \file add.cpp
 * \brief Unit tests for polynomial addition.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-29
 */


# include <random>

# include <gtest/gtest.h>

# include "basic.h"


extern std::default_random_engine generator;

using namespace simpoly;


TEST(PolynomialAdd, EqualLength)
{
    std::uniform_int_distribution<int> dist1(1, 20);
    std::uniform_real_distribution<double> dist2(-5.0, 5.0);
    
    int len = dist1(generator);
    
    std::valarray<double> p1(len), p2(len), result(len);
    
    for(unsigned i=0; i<len; ++i)
    {
        p1[i] = dist2(generator);
        p2[i] = dist2(generator);
    }
    
    result = simpoly::basic::add(p1, p2);
    
    for(unsigned i=0; i<len; ++i)
        ASSERT_NEAR(p1[i]+p2[i], result[i], 1e-12);
}

TEST(PolynomialAdd, FirstOneLonger)
{
    std::uniform_int_distribution<int> dist1(1, 20);
    
    int len = dist1(generator);
    std::valarray<double> p1(len), result(len);
    
    std::uniform_int_distribution<int> dist2(0, len-1);
    std::valarray<double> p2(dist2(generator));
    
    std::uniform_real_distribution<double> rdouble(-5.0, 5.0);
    for(unsigned i=0; i<len; ++i) p1[i] = rdouble(generator);
    for(unsigned i=0; i<p2.size(); ++i) p2[i] = rdouble(generator);
    
    result = simpoly::basic::add(p1, p2); 
    for(unsigned i=0; i<p2.size(); ++i) ASSERT_NEAR(p1[i]+p2[i], result[i], 1e-12);
    for(unsigned i=p2.size(); i<len; ++i) ASSERT_NEAR(p1[i], result[i], 1e-12);
}

TEST(PolynomialAdd, SecondOneLonger)
{
    std::uniform_int_distribution<int> dist1(1, 20);
    
    int len = dist1(generator);
    std::valarray<double> p1(len), result(len);
    
    std::uniform_int_distribution<int> dist2(0, len-1);
    std::valarray<double> p2(dist2(generator));
    
    std::uniform_real_distribution<double> rdouble(-5.0, 5.0);
    for(unsigned i=0; i<len; ++i) p1[i] = rdouble(generator);
    for(unsigned i=0; i<p2.size(); ++i) p2[i] = rdouble(generator);
    
    result = simpoly::basic::add(p2, p1); 
    for(unsigned i=0; i<p2.size(); ++i) ASSERT_NEAR(p1[i]+p2[i], result[i], 1e-12);
    for(unsigned i=p2.size(); i<len; ++i) ASSERT_NEAR(p1[i], result[i], 1e-12);
}

TEST(PolynomialGCD, GCD1)
{
    basic::DArry p1({
        6.48000000e+02,   5.40000000e+02,  -8.10000000e+02,
        -1.11900000e+03,  -1.48000000e+02,   4.43000000e+02,
        3.26000000e+02,   1.03000000e+02,   1.60000000e+01,
        1.00000000e+00});
    
    basic::DArry q1_expect({-108., -108., 45., 104., 54., 12., 1.});
    basic::DArry q2_expect({18., 21., 8., 1.});
    basic::DArry q3_expect({3., 1.});
    basic::DArry q4_expect({1.});
    
    basic::DArry q1 = basic::GCD(p1, basic::derivative(p1));
    basic::DArry q2 = basic::GCD(q1, basic::derivative(q1));
    basic::DArry q3 = basic::GCD(q2, basic::derivative(q2));
    basic::DArry q4 = basic::GCD(q3, basic::derivative(q3));
    
    for(unsigned i=0; i<q1_expect.size(); ++i)
        ASSERT_NEAR(q1_expect[i], q1[i], 1e-8);
    
    for(unsigned i=0; i<q2_expect.size(); ++i)
        ASSERT_NEAR(q2_expect[i], q2[i], 1e-8);
    
    for(unsigned i=0; i<q3_expect.size(); ++i)
        ASSERT_NEAR(q3_expect[i], q3[i], 1e-8);
    
    for(unsigned i=0; i<q4_expect.size(); ++i)
        ASSERT_NEAR(q4_expect[i], q4[i], 1e-8);
    
}

TEST(PolynomialGCD, GCD2)
{
    basic::DArry p1({
            -7.5202173502138876682e-04,  -8.6835988738648978158e-03,
            5.5348919281959091387e-02,   1.4315233176654140745e-01,
            -4.3769681353357492437e-01,  -1.2069048603868517411e+00,
            -2.8638873139566145554e-02,   1.6438835517368903805e+00,
            1.0000000000000000000e+00});
    
    basic::DArry q1_expect({-0.2431327906295455976, 1.0});
    basic::DArry q2_expect({1.});
    
    basic::DArry q1 = basic::GCD(p1, basic::derivative(p1));
    basic::DArry q2 = basic::GCD(q1, basic::derivative(q1));
    
    for(unsigned i=0; i<q1_expect.size(); ++i)
        ASSERT_NEAR(q1_expect[i], q1[i], 1e-8);
    
    for(unsigned i=0; i<q2_expect.size(); ++i)
        ASSERT_NEAR(q2_expect[i], q2[i], 1e-8);
    
}

TEST(PolynomialGCD, GCD3)
{
    basic::DArry p1({
            0.0000000000000000000e+00,  -2.9491200000000012543e-02,
            8.8473600000000013344e-02,   9.6112640000000015839e-01,
            -3.1193088000000002147e+00,  -5.5872000000000010544e+00,
            2.5158400000000003871e+01,  -2.2080000000000006288e+00,
            -6.3264000000000010004e+01,   5.2000000000000014211e+01,
            3.5999999999999992895e+01,  -6.0000000000000000000e+01,
            2.0000000000000000000e+01});
    
    basic::DArry q1_expect({1.0, -2.0, 1.0});
    basic::DArry q2_expect({-1.0, 1.0});
    basic::DArry q3_expect({1.});
    
    basic::DArry q1 = basic::GCD(p1, basic::derivative(p1));
    basic::DArry q2 = basic::GCD(q1, basic::derivative(q1));
    basic::DArry q3 = basic::GCD(q2, basic::derivative(q2));
    
    for(unsigned i=0; i<q1_expect.size(); ++i)
        ASSERT_NEAR(q1_expect[i], q1[i], 1e-8);
    
    for(unsigned i=0; i<q2_expect.size(); ++i)
        ASSERT_NEAR(q2_expect[i], q2[i], 1e-8);
    
    for(unsigned i=0; i<q3_expect.size(); ++i)
        ASSERT_NEAR(q3_expect[i], q3[i], 1e-8);
    
}
