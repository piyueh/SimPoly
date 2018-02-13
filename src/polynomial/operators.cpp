/**
 * \file operators.cpp
 * \brief Implementation of the operators of class Polynomial.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-28
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

// Polynomial::divide
Polynomial Polynomial::divide(const Polynomial &divisor, Polynomial &R)
{
    Polynomial Q;
    Q = poly::divide(*this, divisor, R);
    return std::move(Q);
}

// Polynomial::quotient
Polynomial Polynomial::quotient(const Polynomial &divisor)
{ return std::move(poly::quotient(*this, divisor)); }

// Polynomial::remainder
Polynomial Polynomial::remainder(const Polynomial &divisor)
{ return std::move(poly::remainder(*this, divisor)); }

// evaluation
double Polynomial::operator()(const double x) const
{
    if (_use_roots) return evaluate_from_root(_coef.back(), _rroots, x);
    return evaluate(_coef, x);
}

// evaluation
DArry Polynomial::operator()(const DArry &x) const
{
    DArry result(x.size());
    if (_use_roots)
        for(unsigned i=0; i<x.size(); ++i)
            result[i] = evaluate_from_root(_coef.back(), _rroots, x[i]);

    for(unsigned i=0; i<x.size(); ++i)
        result[i] = evaluate(_coef, x[i]);

    return result;
}

// copy assignment
Polynomial &Polynomial::operator=(const Polynomial &p) = default;

// move assigment
Polynomial &Polynomial::operator=(Polynomial &&p) = default;

// +=
Polynomial & Polynomial::operator+=(const Polynomial &rhs)
{
    set(add(_coef, rhs._coef));
    return *this;
}

// +=
Polynomial & Polynomial::operator+=(const double &rhs)
{
    set(add(_coef, rhs));
    return *this;
}

// -=
Polynomial & Polynomial::operator-=(const Polynomial &rhs)
{
    set(substract(_coef, rhs._coef));
    return *this;
}

// -=
Polynomial & Polynomial::operator-=(const double &rhs)
{
    set(substract(_coef, rhs));
    return *this;
}

// *=
Polynomial & Polynomial::operator*=(const Polynomial &rhs)
{
    set(multiply(_coef, rhs._coef));
    return *this;
}

// *=
Polynomial & Polynomial::operator*=(const double &rhs)
{
    set(multiply(_coef, rhs));
    return *this;
}

// /=
Polynomial & Polynomial::operator/=(const double &rhs)
{
    set(basic::divide(_coef, rhs));
    return *this;
}

// ==
bool Polynomial::operator==(const Polynomial &rhs)
{
    if (_d != rhs._d) return false;

    for(unsigned i=0; i<_coef.size(); ++i)
        if (std::abs(_coef[i] - rhs._coef[i]) > 1e-12) return false;

    return true;
}

// !=
bool Polynomial::operator!=(const Polynomial &rhs)
{
    return ! (this->operator==(rhs));
}

// +
Polynomial operator+(Polynomial lhs, const Polynomial &rhs)
{ return std::move(lhs += rhs); }

// +
Polynomial operator+(Polynomial lhs, const double &rhs)
{ return std::move(lhs += rhs); }

// +
Polynomial operator+(const double &lhs, Polynomial rhs)
{ return std::move(rhs += lhs); }

// -
Polynomial operator-(Polynomial lhs, const Polynomial &rhs)
{ return std::move(lhs -= rhs); }

// -
Polynomial operator-(Polynomial lhs, const double &rhs)
{ return std::move(lhs -= rhs); }

// -
Polynomial operator-(const double &lhs, Polynomial rhs)
{ return std::move(Polynomial(substract(lhs, rhs._coef))); }

// *
Polynomial operator*(Polynomial lhs, const Polynomial &rhs)
{ return std::move(lhs *= rhs); }

// *
Polynomial operator*(Polynomial lhs, const double &rhs)
{ return std::move(lhs *= rhs); }

// *
Polynomial operator*(const double &lhs, Polynomial rhs)
{ return std::move(rhs *= lhs); }

// /
Polynomial operator/(Polynomial lhs, const double &rhs)
{ return std::move(lhs /= rhs); }

// %
Polynomial operator%(Polynomial lhs, const Polynomial &rhs)
{ return std::move(remainder(lhs, rhs)); }

// poly::divide
Polynomial divide(const Polynomial &p1, const Polynomial &p2, Polynomial &R)
{
    DArry q, r;
    q = basic::divide(p1._coef, p2._coef, r);
    R = Polynomial(r);
    return std::move(Polynomial(q));
}

// poly::quotient
Polynomial quotient(const Polynomial &p1, const Polynomial &p2)
{
    DArry q, r;
    q = basic::divide(p1._coef, p2._coef, r);
    return std::move(Polynomial(q));
}

// poly::remainder
Polynomial remainder(const Polynomial &p1, const Polynomial &p2)
{
    DArry q, r;
    q = basic::divide(p1._coef, p2._coef, r);
    return std::move(Polynomial(r));
}

} // end of namespace poly
} // end of namespace simpoly
