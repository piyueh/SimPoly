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
