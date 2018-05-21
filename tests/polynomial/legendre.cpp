/**
 * \file tests/polynomial/legendre.cpp
 * \brief Unit tests for functions creating Legendre polynomials.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-05-21
 */

# include <gtest/gtest.h>

# include "exceptions.h"
# include "polynomial.h"

using namespace simpoly;

TEST(LegendrePoly, LegendrePoly)
{
    std::vector<poly::Polynomial> expect(11);

    expect[0] = poly::Polynomial({1.0});
    expect[1] = poly::Polynomial({0.0, 1.0});
    expect[2] = poly::Polynomial({-0.5, 0.0, 1.5});
    expect[3] = poly::Polynomial({0. , -1.5,  0. ,  2.5});
    expect[4] = poly::Polynomial({0.375, 0.0, -3.75, 0.0, 4.375});
    expect[5] = poly::Polynomial({0.0, 1.875, 0.0, -8.75, 0.0, 7.875});
    expect[6] = poly::Polynomial({-0.3125, 0.0, 6.5625, 0.0, -19.6875, 0.0,
            14.4375});
    expect[7] = poly::Polynomial({0.0, -2.1875, 0.0, 19.6875, 0.0, -43.3125,
            0.0, 26.8125});
    expect[8] = poly::Polynomial({0.2734375, 0.0, -9.84375, 0.0, 54.140625,
            0.0, -93.84375, 0.0, 50.2734375});
    expect[9] = poly::Polynomial({0.0, 2.4609375, 0.0, -36.09375, 0.0,
            140.765625, 0.0, -201.09375, 0.0, 94.9609375});
    expect[10] = poly::Polynomial({-0.24609375, 0.0, 13.53515625, 0.0,
            -117.3046875, 0.0, 351.9140625, 0.0, -427.32421875, 0.0, 180.42578125});

    std::vector<basic::DArry> expect_roots(11);

    expect_roots[0] = {};
    expect_roots[1] = {0.0};
    expect_roots[2] = {0.577350269189626, -0.577350269189626};
    expect_roots[3] = {0.774596669241483, -0.774596669241484,  0.};
    expect_roots[4] = {0.861136311594053, -0.861136311594053, 0.339981043584857,
        -0.339981043584856};
    expect_roots[5] = {-0.906179845938664, -0.538469310105683, 0.906179845938664,
        0.538469310105683, 0.};
    expect_roots[6] = {-0.932469514203151, -0.661209386466264, 0.932469514203152,
        0.661209386466263, -0.238619186083197, 0.238619186083197};
    expect_roots[7] = {0.94910791234276, 0.741531185599395, 0.405845151377397,
        -0.949107912342758, -0.741531185599396, -0.405845151377397, 0.};
    expect_roots[8] = {-0.960289856497528, -0.796666477413634, 0.960289856497534,
        0.796666477413628, -0.525532409916327, 0.52553240991633, -0.18343464249565,
        0.18343464249565};
    expect_roots[9] = {-0.968160239507631, -0.836031107326625, -0.613371432700594,
        -0.324253423403809, 0.968160239507627, 0.836031107326635, 0.613371432700591,
        0.324253423403809, 0.};
    expect_roots[10] = {-0.973906528517121, -0.865063366689068, -0.679409568298991,
        0.973906528517163, 0.865063366688997, 0.67940956829902, -0.433395394129251,
        0.433395394129247, -0.148874338981631, 0.148874338981631};

    for(unsigned degree=0; degree<11; ++degree)
    {
        poly::Polynomial p = poly::Legendre(degree);

        ASSERT_EQ(expect[degree], p);
        ASSERT_EQ(poly::PolyType::LEGENDRE, p.type());
        ASSERT_NEAR(1.0, p(1.0), 1e-12);
        ASSERT_NEAR(((degree%2==0)?1.0:-1.0), p(-1.0), 1e-12);

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
