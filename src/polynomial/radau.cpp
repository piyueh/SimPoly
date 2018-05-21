/**
 * \file radau.cpp
 * \brief Implementation of factory creating Radau polynomials.
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

// function to create Radau polynomials
Polynomial Radau(const unsigned n, const PolyType type)
{
    Polynomial s;
    Polynomial p;
    double scaling;

    // end node polynomial
    switch (type)
    {
        case PolyType::RIGHTRADAU:
            s = Jacobi(1, 0, n-1);
            p = Polynomial({1.0, -1.0});
            scaling = 0.5 * ((n%2==0)?-1.0:1.0);
            break;
        case PolyType::LEFTRADAU:
            s = Jacobi(0, 1, n-1);
            p = Polynomial({1.0, 1.0});
            scaling = 0.5;
            break;
        default:
            throw exceptions::IllegalType(__FL__, "RIGHTRADAU and LEFTRADAU");
    }

    p *= s;
    p *= scaling;
    p.set(type);

    return std::move(p);
}

} // end of namespace poly
} // end of namespace simpoly
