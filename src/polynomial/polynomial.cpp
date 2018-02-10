/**
 * \file polynomial.cpp
 * \brief Implementation of the class Polynomial.
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

// copy constructor
Polynomial::Polynomial(const Polynomial &p) = default;

// move constructor
Polynomial::Polynomial(Polynomial &&p) = default;

// constructor
Polynomial::Polynomial(const DArry &coef, const PolyType &type):
    _type(type) { set(coef); }

// constructor
Polynomial::Polynomial(const double l, const DArry &roots, const PolyType &type):
    _type(type) { set(l, roots); }

// constructor
Polynomial::Polynomial(const double l, const CArry &roots, const PolyType &type):
    _type(type) { set(l, roots); }

// constructor
Polynomial::Polynomial(const double l, const DArry &rroots, const CArry &croots,
        const PolyType &type): _type(type) { set(l, rroots, croots); }

// constructor
Polynomial::Polynomial(const DArry &coef, const DArry &roots,
        const PolyType &type): _type(type) { set(coef, roots); }

// constructor
Polynomial::Polynomial(const DArry &coef, const CArry &roots,
        const PolyType &type): _type(type) { set(coef, roots); }

// constructor
Polynomial::Polynomial(const DArry &coef, const DArry &rroots, const CArry &croots,
        const PolyType &type): _type(type) { set(coef, rroots, croots); }


// re-set type
void Polynomial::set(const PolyType &type) { _type = type; }

// re-set coefficient
void Polynomial::set(const DArry &coef, const double &tol)
{
    CArry tmp;

    _coef = coef;// copy/move coef to _coef
    _d = _coef.size() - 1; // get degree of polynomial
    _nrr = _ncr = 0; // initialize _nrr and _ncr
    _rroots.clear();
    _croots.clear();
    _use_roots = false;

    tmp = yan_and_chieng_2006(_coef);// get roots
    for(const auto &it: tmp)
    {
        if (std::abs(it.imag()) < tol)
            _rroots.push_back(it.real());
        else
            _croots.push_back(it);
    }

    _nrr = _rroots.size();
    _ncr = _croots.size();
}

// re-set roots
void Polynomial::set(const double l, const DArry &roots)
{
    _rroots = roots;
    _d = _nrr = _rroots.size();
    _use_roots = true;

    _croots.clear();
    _ncr = 0;

    _coef = to_coefficients(l, _rroots);
}

// re-set roots
void Polynomial::set(const double l, const CArry &roots)
{
    _croots = roots;
    _d = _ncr = _croots.size();
    _use_roots = false;

    _rroots.clear();
    _nrr = 0;

    // note this Polynomial class can only handle real-number coefficients
    _coef = to_DArry(to_coefficients(Cmplx(l), _croots));
}

// re-set roots
void Polynomial::set(const double l, const DArry &rroots, const CArry &croots)
{
    _rroots = rroots;
    _croots = croots;

    _nrr = _rroots.size();
    _ncr = _croots.size();
    _d = _nrr + _ncr;

    _use_roots = false;

    // note this Polynomial class can only handle real-number coefficients
    _coef = to_DArry(to_coefficients(Cmplx(l), _croots));
    _coef = multiply(_coef, to_coefficients(1.0, _rroots));
}

// re-set both roots and coefficients
void Polynomial::set(const DArry &coef, const DArry &roots)
{
    _coef = coef;
    _d = _coef.size() - 1;

    _rroots = roots;
    _nrr = _rroots.size();

    _croots.clear();
    _ncr = 0;

    _use_roots = true;

# ifndef NDEBUG
    if (_d != _nrr) throw exceptions::UnmatchedLength(__FL__, _d, _nrr);

    for(const auto &it: _rroots)
    {
        double value = evaluate(_coef, it);
        if (std::abs(value) > 1e-12) throw exceptions::ExpectingZero(__FL__, value);
    }
# endif
}

// re-set both roots and coefficients
void Polynomial::set(const DArry &coef, const CArry &roots)
{
    _coef = coef;
    _d = _coef.size() - 1;

    _croots = roots;
    _ncr = _croots.size();

    _rroots.clear();
    _nrr = 0;

    _use_roots = false;

# ifndef NDEBUG
    if (_d != _ncr) throw UnmatchedLength(__FL__, _d, _ncr);

    for(const auto &it: _croots)
    {
        Cmplx value = evaluate(to_CArry(_coef), it);
        if (std::abs(value.imag()) > 1e-12) throw FoundComplexNumber(__FL__, value);
        if (std::abs(value.real()) > 1e-12) throw ExpectingZero(__FL__, value.real());
    }
# endif
}

// re-set both roots and coefficients
void Polynomial::set(const DArry &coef,
        const DArry &rroots, const CArry &croots)
{
    _coef = coef;
    _rroots = rroots;
    _croots = croots;

    _d = _coef.size() - 1;
    _nrr = _rroots.size();
    _ncr = _croots.size();

    _use_roots = false;

# ifndef NDEBUG
    if (_d != (_nrr + _ncr)) throw UnmatchedLength(__FL__, _d, _ncr+_nrr);

    for(const auto &it: _rroots)
    {
        double value = evaluate(_coef, it);
        if (std::abs(value) > 1e-12) throw ExpectingZero(__FL__, value);
    }

    for(const auto &it: _croots)
    {
        Cmplx value = evaluate(to_CArry(_coef), it);
        if (std::abs(value.imag()) > 1e-12) throw FoundComplexNumber(__FL__, value);
        if (std::abs(value.real()) > 1e-12) throw ExpectingZero(__FL__, value.real());
    }
# endif
}

// re-set only one coefficient
void Polynomial::set(const int d, const double value)
{
    _coef[d] = value;
    set(_coef);
}

// get a string of type
PolyType Polynomial::type() const { return _type; }

// get a bool indicating whether to use roots or coefficients for evaluation
bool Polynomial::use_roots() const { return _use_roots; }

// get degree of polynomial
unsigned Polynomial::degree() const { return _d; }

// get number of real roots
unsigned Polynomial::n_real_roots() const { return _nrr; }

// get number of complex roots
unsigned Polynomial::c_real_roots() const { return _ncr; }

// get coefficients
DArry Polynomial::coef() const { return _coef; }

// get real roots
DArry Polynomial::real_roots() const { return _rroots; }

// get complex roots
CArry Polynomial::cmplx_roots() const { return _croots; }

// get all roots in a complex vector
CArry Polynomial::roots() const
{
    CArry result(_d);
    std::copy(_rroots.begin(), _rroots.end(), result.begin());
    std::copy(_croots.begin(), _croots.end(), result.begin()+_nrr);
    return result;
}

// derivative
Polynomial Polynomial::deriv() const { return Polynomial(derivative(_coef)); }

// integral
Polynomial Polynomial::integ() const { return Polynomial(integral(_coef)); }

// evaluation
double Polynomial::operator()(const double x) const
{
    if (_use_roots) return evaluate_from_root(_coef.back(), _rroots, x);
    return evaluate(_coef, x);
}

// copy assignment
Polynomial &Polynomial::operator=(const Polynomial &p) = default;

// move assigment
Polynomial &Polynomial::operator=(Polynomial &&p) = default;

} // end of namespace poly
} // end of namespace simpoly
