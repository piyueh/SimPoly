/**
 * \file tests/polynomial/radau.cpp
 * \brief Unit tests for functions creating Radau polynomials.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-05-21
 */

# include <gtest/gtest.h>

# include "exceptions.h"
# include "polynomial.h"

using namespace simpoly;

// TODO: what happen if users pass n=0 to Radau function?

TEST(RadauPoly, LeftRadauPoly)
{
    std::vector<poly::Polynomial> expect(11);

    expect[0] = poly::Polynomial(); // degree 1 Radau doesn't exist
    expect[1] = poly::Polynomial({0.5, 0.5});
    expect[2] = poly::Polynomial({-0.25, 0.5, 0.75});
    expect[3] = poly::Polynomial({-0.25, -0.75, 0.75, 1.25});
    expect[4] = poly::Polynomial({0.1875, -0.75, -1.875, 1.25, 2.1875});
    expect[5] = poly::Polynomial({0.1875, 0.9375, -1.875, -4.375, 2.1875, 3.9375});
    expect[6] = poly::Polynomial({
            -0.15625, 0.9375, 3.28125, -4.375, -9.84375, 3.9375, 7.21875});
    expect[7] = poly::Polynomial({
            -0.15625, -1.09375, 3.28125, 9.84375,
            -9.84375, -21.65625, 7.21875, 13.40625 });
    expect[8] = poly::Polynomial({
            0.13671875, -1.09375, -4.921875, 9.84375, 27.0703125,
            -21.65625, -46.921875, 13.40625, 25.13671875});
    expect[9] = poly::Polynomial({
            0.13671875, 1.23046875, -4.921875, -18.046875, 27.0703125,
            70.3828125, -46.921875, -100.546875, 25.13671875, 47.48046875});
    expect[10] = poly::Polynomial({
            -1.23046875e-01, 1.23046875e+00, 6.767578125e+00, -1.8046875e+01,
            -5.865234375e+01, 7.03828125e+01, 1.7595703125e+02, -1.00546875e+02,
            -2.13662109375e+02, 4.748046875e+01, 9.0212890625e+01});

    std::vector<basic::DArry> expect_roots(11);

    expect_roots[0] = {}; // degree 1 Radau doesn't exist
    expect_roots[1] = {-1.};
    expect_roots[2] = {-1., 0.333333333333333};
    expect_roots[3] = {0.689897948556635, -1., -0.289897948556636};
    expect_roots[4] = {
        0.822824080974592, -1., -0.575318923521693, 0.181066271118531};
    expect_roots[5] = {
        0.885791607770964, -1., -0.720480271312438, 0.446313972723753,
        -0.167180864737834};
    expect_roots[6] = {
        0.920380285897062, -1., -0.802929828402348, 0.603973164252784,
        -0.390928546707273, 0.124050379505228};
    expect_roots[7] = {
        0.941367145680432, 0.703842800663031, -1., -0.853891342639482,
        -0.538467724060108, 0.326030619437691, -0.1173430375431};
    expect_roots[8] = {
        0.955041227122574, 0.770641893678193, -1., -0.887474878926158,
        -0.639518616526216, 0.468420354430821, -0.294750565773661,
        0.094307252661111};
    expect_roots[9] = {
        0.964440169705265, 0.817352784200416, 0.571383041208737, -1.,
        -0.910732089420092, -0.711267485915699, -0.426350485711141,
        0.256135670833456, -0.090373369606853};
    expect_roots[10] = {
        0.971175180702248, 0.851225220581606, -1., -0.927484374233556,
        -0.763842042419999, 0.647766687674011, -0.525646030370082,
        0.380664840144724, -0.236234469390588, 0.076059197837978};

    // start from degree 1. 0-degree Radau doesn't make sense
    for(unsigned degree=1; degree<11; ++degree)
    {
        poly::Polynomial p = poly::Radau(degree, poly::PolyType::LEFTRADAU);

        ASSERT_EQ(expect[degree], p)<< std::setprecision(15) ;
        ASSERT_EQ(poly::PolyType::LEFTRADAU, p.type());
        ASSERT_NEAR(0.0, p(-1.0), 1e-12);
        ASSERT_NEAR(1.0, p(1.0), 1e-12);
        ASSERT_NEAR(((degree%2==0)?-1.0:1.0)*double(degree)/2.0, p.deriv()(-1.0), 1e-12);
        ASSERT_NEAR(double(degree)*double(degree)/2.0, p.deriv()(1.0), 1e-12);

        basic::DArry r_roots = p.real_roots();
        basic::CArry c_roots = p.cmplx_roots();
        unsigned nr = p.n_real_roots();
        unsigned nc = p.n_cmplx_roots();

        ASSERT_EQ(degree, r_roots.size());
        ASSERT_EQ(degree, nr);
        ASSERT_EQ(0, c_roots.size());
        ASSERT_EQ(0, nc);

        std::sort(std::begin(expect_roots[degree]), std::end(expect_roots[degree]));
        std::sort(std::begin(r_roots), std::end(r_roots));

        for(unsigned i=0; i<nr; ++i)
            ASSERT_NEAR(expect_roots[degree][i], r_roots[i], 1e-10);
    }
}

TEST(RadauPoly, RightRadauPoly)
{
    std::vector<poly::Polynomial> expect(11);

    expect[0] = poly::Polynomial(); // degree 1 Radau doesn't exist
    expect[1] = poly::Polynomial({
            0.5, -0.5});
    expect[2] = poly::Polynomial({
            -0.25, -0.5, 0.75});
    expect[3] = poly::Polynomial({
            -0.25, 0.75, 0.75, -1.25});
    expect[4] = poly::Polynomial({
            0.1875, 0.75, -1.875, -1.25, 2.1875});
    expect[5] = poly::Polynomial({
            0.1875, -0.9375, -1.875, 4.375, 2.1875, -3.9375});
    expect[6] = poly::Polynomial({
            -0.15625, -0.9375, 3.28125, 4.375, -9.84375, -3.9375, 7.21875});
    expect[7] = poly::Polynomial({
            -0.15625, 1.09375, 3.28125, -9.84375, -9.84375, 21.65625, 7.21875,
            -13.40625});
    expect[8] = poly::Polynomial({
            0.13671875, 1.09375, -4.921875, -9.84375, 27.0703125, 21.65625,
            -46.921875, -13.40625, 25.13671875});
    expect[9] = poly::Polynomial({
            0.13671875, -1.23046875, -4.921875, 18.046875, 27.0703125,
            -70.3828125, -46.921875, 100.546875, 25.13671875, -47.48046875});
    expect[10] = poly::Polynomial({
            -1.23046875e-01, -1.23046875e+00, 6.767578125e+00, 1.8046875e+01,
            -5.865234375e+01, -7.03828125e+01, 1.7595703125e+02, 1.00546875e+02,
            -2.13662109375e+02, -4.748046875e+01, 9.0212890625e+01});

    std::vector<basic::DArry> expect_roots(11);

    expect_roots[0] = {}; // degree 1 Radau doesn't exist
    expect_roots[1] = {
        1.};
    expect_roots[2] = {
        1., -0.333333333333333};
    expect_roots[3] = {
        -0.689897948556636, 1., 0.289897948556636};
    expect_roots[4] = {
        -0.822824080974593, 1., 0.575318923521694, -0.181066271118531};
    expect_roots[5] = {
        -0.885791607770965, 1., 0.720480271312439, -0.446313972723752,
        0.167180864737834};
    expect_roots[6] = {
        -0.920380285897063, 1., 0.802929828402348, -0.603973164252783,
        0.390928546707272, -0.124050379505228};
    expect_roots[7] = {
        -0.941367145680432, -0.703842800663031, 1., 0.853891342639486,
        0.538467724060109, -0.326030619437691, 0.1173430375431};
    expect_roots[8] = {
        -0.955041227122572, -0.770641893678195, 1., 0.887474878926159,
        0.639518616526215, -0.46842035443082, 0.294750565773661,
        -0.094307252661111};
    expect_roots[9] = {
        -0.964440169705274, -0.817352784200411, -0.571383041208739, 1.,
        0.910732089420068, 0.711267485915706, 0.42635048571114,
        -0.256135670833455, 0.090373369606853};
    expect_roots[10] = {
        -0.971175180702222, -0.851225220581646, 1., 0.927484374233571,
        0.763842042420008, -0.647766687673994, 0.525646030370078,
        -0.380664840144725, 0.236234469390588, -0.076059197837978};


    // start from degree 1. 0-degree Radau doesn't make sense
    for(unsigned degree=1; degree<11; ++degree)
    {
        poly::Polynomial p = poly::Radau(degree, poly::PolyType::RIGHTRADAU);

        ASSERT_EQ(expect[degree], p)<< std::setprecision(15) ;
        ASSERT_EQ(poly::PolyType::RIGHTRADAU, p.type());
        ASSERT_NEAR(1.0, p(-1.0), 1e-12);
        ASSERT_NEAR(0.0, p(1.0), 1e-12);
        ASSERT_NEAR(-double(degree)*double(degree)/2.0, p.deriv()(-1.0), 1e-12);
        ASSERT_NEAR(((degree%2==0)?1.0:-1.0)*double(degree)/2.0, p.deriv()(1.0), 1e-12);

        basic::DArry r_roots = p.real_roots();
        basic::CArry c_roots = p.cmplx_roots();
        unsigned nr = p.n_real_roots();
        unsigned nc = p.n_cmplx_roots();

        ASSERT_EQ(degree, r_roots.size());
        ASSERT_EQ(degree, nr);
        ASSERT_EQ(0, c_roots.size());
        ASSERT_EQ(0, nc);

        std::sort(std::begin(expect_roots[degree]), std::end(expect_roots[degree]));
        std::sort(std::begin(r_roots), std::end(r_roots));

        for(unsigned i=0; i<nr; ++i)
            ASSERT_NEAR(expect_roots[degree][i], r_roots[i], 1e-10);
    }
}
