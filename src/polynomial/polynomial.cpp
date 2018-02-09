/**
 * \file polynomial.cpp
 * \brief Implementation of the class Polynomial.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-28
 */


# include "polynomial.h"


using namespace simpoly::basic;


# define NotImplemented \
    throw std::string("Function in ") + __FILE__ + \
        "(line " + std::to_string(__LINE__) +") not implemented yet.";

namespace simpoly
{
namespace poly
{
    
Polynomial::Polynomial(Polynomial &&p) { NotImplemented; }

Polynomial::Polynomial(const DArry &coef): _coef(coef) { NotImplemented; }

Polynomial::Polynomial(const double l, const DArry &rroots) { NotImplemented; }

Polynomial::Polynomial(const double l, const CArry &croots) { NotImplemented; }

Polynomial::Polynomial(const double l, 
        const DArry &rroots, const CArry &croots) { NotImplemented; }

Polynomial::Polynomial(const DArry &coef, const DArry &rroots) { NotImplemented; }

Polynomial::Polynomial(const DArry &coef, const CArry &croots) { NotImplemented; }

Polynomial::Polynomial(const DArry &coef, 
        const DArry &roots, const CArry &croots) { NotImplemented; }


void Polynomial::set(const DArry &coef) { NotImplemented; }

void Polynomial::set(const double l, const DArry &r) { NotImplemented; }

void Polynomial::set(const double l, const CArry &r) { NotImplemented; }

void Polynomial::set(const double l, const DArry &rr, const CArry &rc) { NotImplemented; }

void Polynomial::set(const DArry &coef, const DArry &rroots) { NotImplemented; }

void Polynomial::set(const DArry &coef, const CArry &croots) { NotImplemented; }

void Polynomial::set(const DArry &coef, 
        const DArry &roots, const CArry &croots) { NotImplemented; }

void Polynomial::set(const int d, const double value) { NotImplemented; }


unsigned Polynomial::degree() const { NotImplemented; }

DArry Polynomial::coef() const { NotImplemented; }

CArry Polynomial::roots() const { NotImplemented; }

DArry Polynomial::real_roots() const { NotImplemented; }

CArry Polynomial::cmplx_roots() const { NotImplemented; }


Polynomial Polynomial::deriv() const { NotImplemented; }

Polynomial Polynomial::integ() const { NotImplemented; }


double Polynomial::operator()(const double x) const { NotImplemented; }

    
} // end of namespace poly
} // end of namespace simpoly
