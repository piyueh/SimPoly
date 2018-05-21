/**
 * \file jacobi.cpp
 * \brief Implementation of factory creating Jacobi-family polynomials.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-02-14
 */

# include <cmath>

# include "exceptions.h"
# include "polynomial.h"


using namespace simpoly::basic;
using namespace simpoly::exceptions;


namespace simpoly
{
namespace poly
{

// function to create Jacobi-family polynomials
Polynomial Jacobi(const double alpha, const double beta, const unsigned n)
{
# ifndef NDEBUG
    if (alpha <= -1.0) throw JacobiParameters(__FL__, alpha, beta);
    if (beta <= -1.0) throw JacobiParameters(__FL__, alpha, beta);
# endif

    // P_{0}
    if (n == 0) return Polynomial({1.0});

    // P_{1}
    if (n == 1) return Polynomial({(alpha-beta)/2., (alpha+beta)/2.+1.});

    // previous and current polynomial
    Polynomial Pi({(alpha-beta)/2., (alpha+beta)/2.+1.}); // P_i = P_1
    Polynomial Pim1({1.0}); // P_{i-1} == P_0

    const double &c1 = alpha + beta; // for convience

    // the result of each loop is P_{i+1}
    for(unsigned i=1; i<n; ++i)
    {
        double np1 = i + 1;
        double nt2 = i * 2;
        double np1t2 = (i + 1) * 2;
        double nt2p1 = i * 2 + 1;

        double a1 = np1t2 * (np1 + c1) * (nt2 + c1);
        double a2 = (nt2p1 + c1) * c1 * (alpha - beta);
        double a3 = (nt2 + c1) * (nt2p1 + c1) * (np1t2 + c1);
        double a4 = 2.0 * (alpha + i) * (beta + i) * (c1 + np1t2);

        // we re-use the memory space of Pim1 to store P_{i+1}
        Pim1 *= (-a4);
        Pim1 += (Polynomial({a2, a3}) * Pi);
        Pim1 /= a1; // now Pim1 becomes P_{i+1}

        // swap
        std::swap(Pim1, Pi);
    }

    Pi.set(PolyType::JACOBI);

    return std::move(Pi);
}

} // end of namespace poly
} // end of namespace simpoly
