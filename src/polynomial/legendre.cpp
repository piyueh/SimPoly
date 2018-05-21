/**
 * \file legendre.cpp
 * \brief Implementation of factory creating Legendre polynomials.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-05-21
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

// function to create Legendre polynomials
Polynomial Legendre(const unsigned n)
{
    // for low-degree Legendre, use exact values
    switch (n)
    {
        case 0:
        {
            poly::Polynomial Pi({1.0L});
            Pi.set(PolyType::LEGENDRE);
            return std::move(Pi);
        }
        case 1:
        {
            poly::Polynomial Pi({0.0L, 1.0L});
            Pi.set(PolyType::LEGENDRE);
            return std::move(Pi);
        }
        case 2:
        {
            poly::Polynomial Pi({-0.5L, 0.0L, 1.5L});
            Pi.set(PolyType::LEGENDRE);
            return std::move(Pi);
        }
        case 3:
        {
            poly::Polynomial Pi({0.0L, -1.5L, 0.0L, 2.5L});
            Pi.set(PolyType::LEGENDRE);
            return std::move(Pi);
        }
        case 4:
        {
            poly::Polynomial Pi({0.375L, 0.0L, -3.75L, 0.0L, 4.375L});
            Pi.set(PolyType::LEGENDRE);
            return std::move(Pi);
        }
        case 5:
        {
            poly::Polynomial Pi({0.0L, 1.875L, 0.0L, -8.75L, 0.0L, 7.875L});
            Pi.set(PolyType::LEGENDRE);
            return std::move(Pi);
        }
        case 6:
        {
            poly::Polynomial Pi({-0.3125L, 0.0L, 6.5625L, 0.0L, -19.6875L, 0.0L, 14.4375L});
            Pi.set(PolyType::LEGENDRE);
            return std::move(Pi);
        }
        default: break;
    }

    // previous and current polynomial
    Polynomial Pi = Legendre(6); // P_i = P_6
    Polynomial Pim1 = Legendre(5); // P_{i-1} == P_5

    // the result of each loop is P_{i+1}
    for(unsigned i=6; i<n; ++i)
    {
        long double np1 = i + 1.0L;
        long double nt2 = i * 2.0L;
        long double np1t2 = (i + 1.0L) * 2.0L;
        long double nt2p1 = i * 2.0L + 1.0L;

        long double a1 = np1t2 * np1 * nt2;
        long double a3 = nt2 * nt2p1 * np1t2;
        long double a4 = 2.0 * i * i * np1t2;

        // we re-use the memory space of Pim1 to store P_{i+1}
        Pim1 *= (-a4);
        Pim1 += (Polynomial({0.0, double(a3)}) * Pi);
        Pim1 /= a1; // now Pim1 becomes P_{i+1}

        // swap
        std::swap(Pim1, Pi);
    }

    Pi.set(PolyType::LEGENDRE);

    return std::move(Pi);
}

} // end of namespace poly
} // end of namespace simpoly
