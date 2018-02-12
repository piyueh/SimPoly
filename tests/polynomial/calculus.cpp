/**
 * \file tests/polynomial/calculus.cpp
 * \brief Unit tests for calculus functions of Polynomial class.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-02-08
 */

# include <gtest/gtest.h>

# include "polynomial.h"

using namespace simpoly;

TEST(PolynomialCalculus, Derivative)
{
    poly::Polynomial p({
        -3.0582015497062718776e-05,   1.2877895859209613928e-03,  -6.8407333802847326235e-04,
        -3.1171792808804225233e-02,   1.1661017762089756836e-02,   2.6233713034532663011e-01,
        -3.2802450138374281141e-02,  -8.7659322558118857316e-01,   1.2229768187730361007e-02,
        1.0000000000000000000e+00});

    basic::DArry expect({
        1.2877895859209613928e-03,  -1.3681466760569465247e-03,  -9.3515378426412679169e-02,
        4.6644071048359027343e-02,   1.3116856517266330950e+00,  -1.9681470083024568685e-01,
        -6.1361525790683195680e+00,   9.7838145501842888052e-02,   9.0000000000000000000e+00});

    poly::Polynomial d = p.deriv();

    basic::DArry result = d.coef();

    for(unsigned i=0; i<expect.size(); i++)
        ASSERT_NEAR(expect[i], result[i], 1e-10);
}

TEST(PolynomialCalculus, Integral)
{
    poly::Polynomial p({
        -1.3883393317951405604e-03, -4.2579391363204877052e-02, -3.3604182219819818400e-03,
        6.2437819297829388887e-01, -6.3438730299142587210e-01, -1.6525521770331836890e+00,
        2.5659233014056495037e+00,  5.6875799592633724533e-01, -2.4247744334853669201e+00,
        1.0000000000000000000e+00});

    basic::DArry expect({
        0.                   , -0.0013883393317951406, -0.0212896956816024385,
        -0.0011201394073273273,  0.1560945482445734722, -0.1268774605982851689,
        -0.2754253628388639297,  0.3665604716293784926,  0.0710947494907921557,
        -0.2694193814983740776,  0.1000000000000000056});

    poly::Polynomial d = p.integ();

    basic::DArry result = d.coef();

    for(unsigned i=0; i<expect.size(); i++)
        ASSERT_NEAR(expect[i], result[i], 1e-10);
}
