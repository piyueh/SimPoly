/**
 * \file tests/polynomial/jacobi.cpp
 * \brief Unit tests for functions creating Jacobi-family polynomial.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-02-13
 */

# include <gtest/gtest.h>

# include "exceptions.h"
# include "polynomial.h"

using namespace simpoly;

# ifndef NDEBUG
TEST(JacobiPoly, TestExceptions)
{
    ASSERT_THROW(poly::Jacobi(-1., 0., 0), exceptions::JacobiParameters);
    ASSERT_THROW(poly::Jacobi(0., -1., 0), exceptions::JacobiParameters);
}
# endif

TEST(JacobiPoly, JacobiPoly1)
{
    poly::Polynomial p = poly::Jacobi(3.6397070770505078, 1.0266395951667127, 13);

    poly::Polynomial expect = poly::Polynomial({
        8.6647962638650988243e-01,  -2.8745036096183085306e+00,  -1.0888669504296072432e+02,
        1.8219208319067512036e+01,   2.0746680106383932980e+03,   8.8687766109362178213e+02,
        -1.3721298111564616192e+04,  -9.1200309736587332736e+03,   3.9327879955206204613e+04,
        3.0924979079538272345e+04,  -5.0181920988314173883e+04,  -4.3032838430113879440e+04,
        2.3323274340020372620e+04,   2.1055113354500550486e+04});

    ASSERT_EQ(expect, p);
    ASSERT_EQ(poly::PolyType::JACOBI, p.type());
}

TEST(JacobiPoly, JacobiPoly2)
{
    poly::Polynomial p = poly::Jacobi(3.325420758348053, -0.42502377411185854, 12);

    poly::Polynomial expect = poly::Polynomial({
        -5.7167036695211337349e-01, -4.3233893850837583628e+00, 4.6431398631484185557e+01,
        1.9200228940711082259e+02, -5.1259357256573605355e+02, -2.0175045235751217660e+03,
        1.4591276878910300638e+03, 7.8095424518275585797e+03, 2.3735601487973534063e+02,
        -1.2320647048811981222e+04, -4.8729193574906321373e+03, 6.7055533658185704553e+03,
        4.0080148272714327504e+03});

    ASSERT_EQ(expect, p);
    ASSERT_EQ(poly::PolyType::JACOBI, p.type());
}

TEST(JacobiPoly, JacobiPoly3)
{
    poly::Polynomial p1 = poly::Jacobi(0., 0., 7);
    poly::Polynomial p2 = poly::Jacobi(0., 0., 7);
    poly::Polynomial p3 = (p1*p2).integ();
    ASSERT_NEAR(2.0/(2.0*7+1), p3(1.0)-p3(-1.0), 1e-10);
}

TEST(JacobiPoly, JacobiPoly4)
{
    poly::Polynomial p1 = poly::Jacobi(0., 0., 5);
    poly::Polynomial p2 = poly::Jacobi(0., 0., 18);
    poly::Polynomial p3 = (p1*p2).integ();
    ASSERT_NEAR(0.0, p3(1.0)-p3(-1.0), 1e-10);
}
