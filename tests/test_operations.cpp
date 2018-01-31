/*
 * test_operations.cpp
 * Copyright (C) 2018 Pi-Yueh Chuang <pychuang@gwu.edu>
 *
 * Distributed under terms of the MIT license.
 */


# include <valarray>
# include <random>
# include <chrono>

# include <gtest/gtest.h>

# include "operations.h"
# include "exceptions.h"


TEST(PolynomialAdd, EqualLength)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> dist1(1, 20);
    std::uniform_real_distribution<double> dist2(-5.0, 5.0);
    
    int degree = dist1(generator);
    
    std::valarray<double> p1(degree), p2(degree), result(degree);
    
    for(unsigned i=0; i<degree; ++i)
    {
        p1[i] = dist2(generator);
        p2[i] = dist2(generator);
    }
    
    result = polynomial::op::add(p1, p2);
    
    for(unsigned i=0; i<degree; ++i)
        ASSERT_NEAR(p1[i]+p2[i], result[i], 1e-12);
}


TEST(PolynomialAdd, FirstOneLonger)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> dist1(1, 20);
    
    int degree = dist1(generator);
    std::valarray<double> p1(degree), result(degree);
    
    std::uniform_int_distribution<int> dist2(0, degree-1);
    std::valarray<double> p2(dist2(generator));
    
    std::uniform_real_distribution<double> rdouble(-5.0, 5.0);
    for(unsigned i=0; i<degree; ++i) p1[i] = rdouble(generator);
    for(unsigned i=0; i<p2.size(); ++i) p2[i] = rdouble(generator);
    
    result = polynomial::op::add(p1, p2); 
    for(unsigned i=0; i<p2.size(); ++i) ASSERT_NEAR(p1[i]+p2[i], result[i], 1e-12);
    for(unsigned i=p2.size(); i<degree; ++i) ASSERT_NEAR(p1[i], result[i], 1e-12);
}


TEST(PolynomialAdd, SecondOneLonger)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator(seed);
    std::uniform_int_distribution<int> dist1(1, 20);
    
    int degree = dist1(generator);
    std::valarray<double> p1(degree), result(degree);
    
    std::uniform_int_distribution<int> dist2(0, degree-1);
    std::valarray<double> p2(dist2(generator));
    
    std::uniform_real_distribution<double> rdouble(-5.0, 5.0);
    for(unsigned i=0; i<degree; ++i) p1[i] = rdouble(generator);
    for(unsigned i=0; i<p2.size(); ++i) p2[i] = rdouble(generator);
    
    result = polynomial::op::add(p2, p1); 
    for(unsigned i=0; i<p2.size(); ++i) ASSERT_NEAR(p1[i]+p2[i], result[i], 1e-12);
    for(unsigned i=p2.size(); i<degree; ++i) ASSERT_NEAR(p1[i], result[i], 1e-12);
}


# ifndef NDEBUG
TEST(PolynomialMult, ZeroLengthTests)
{
    std::valarray<double> p1(0), p2(1);
    
    try { polynomial::op::multiply(p1, p2); }
    catch (polynomial::exceptions::ZeroCoeffsLength e) {};
    
    try { polynomial::op::multiply(p2, p1); }
    catch (polynomial::exceptions::ZeroCoeffsLength e) {};
}
# endif


int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
