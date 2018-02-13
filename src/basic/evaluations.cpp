/**
 * \file evaluations.cpp
 * \brief Polynomial evaluation functions.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-29
 */


# include "basic.h"
# include "exceptions.h"


namespace simpoly
{
namespace basic
{

template <typename T>
T evaluate(const T* const &bg, const int len, const T x)
{
# ifndef NDEBUG
    if (len == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
    if (len < 0) throw exceptions::NegativeCoeffsLength(__FILE__, __LINE__, len);
# endif

    if (len == 1)
        return *bg;
    else
        return *bg + x * evaluate(bg+1, len-1, x);
}


template <typename T>
T evaluate(const Arry<T> &coeffs, const T x)
{
    CHECK_COEFS(coeffs, 1e-12);

    return evaluate(&coeffs[0], coeffs.size(), x);
}


template <typename T>
T evaluate_from_root(const T l, const T* const &bg, const int degree, const T x)
{
# ifndef NDEBUG
    if (degree < 0) throw exceptions::NegativeDegree(__FILE__, __LINE__, degree);
# endif

    if (degree == 0)
        return l;
    else
        return (x - *bg) * evaluate_from_root(l, bg+1, degree-1, x);
}


template <typename T>
T evaluate_from_root(const T l, const Arry<T> &roots, const T x)
{
    return evaluate_from_root(l, &roots[0], roots.size(), x);
}


// explicit instantiation
template double evaluate(const double* const &bg, const int len, const double x);
template Cmplx evaluate(const Cmplx* const &bg, const int len, const Cmplx x);
template double evaluate(const DArry &coeffs, const double x);
template Cmplx evaluate(const CArry &coeffs, const Cmplx x);
template double evaluate_from_root(const double l,
        const double* const &bg, const int degree, const double x);
template Cmplx evaluate_from_root(const Cmplx l,
        const Cmplx* const &bg, const int degree, const Cmplx x);
template double evaluate_from_root(
        const double l, const DArry &roots, const double x);
template Cmplx evaluate_from_root(
        const Cmplx l, const CArry &roots, const Cmplx x);

} // end of namespace basic
} // end of namespace simpoly
