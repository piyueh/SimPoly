/**
 * \file operations.cpp
 * \brief Basic operations of polynomials defined by a std::valarray.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-29
 */

# include <algorithm>
# include <numeric>
# include <cmath>

# include "basic.h"
# include "exceptions.h"


namespace simpoly
{
namespace basic
{

template <typename T>
Arry<T> add(const Arry<T> &p1, const Arry<T> &p2)
{
    CHECK_COEFS(p1, 1e-12);
    CHECK_COEFS(p2, 1e-12);

    const Arry<T> *pl = &p1, *ps = &p2;
    Arry<T> result;

    if (p1.size() < p2.size()) std::swap(ps, pl);

    result.resize(pl->size());

    std::transform(ps->begin(), ps->end(), pl->begin(), result.begin(),
            [](const T &x1, const T &x2)->T{return x1+x2;});

    std::copy(pl->begin()+ps->size(), pl->end(), result.begin()+ps->size());

    // eliminate zero leading coefficients
    trim_leading_zeros(result, 1e-12);

    return result;
}

template <typename T>
Arry<T> add(const Arry<T> &p, const T &c)
{
    CHECK_COEFS(p, 1e-12);

    Arry<T> result(p); result[0] += c;
    return result;
}

template <typename T>
Arry<T> add(const T &c, const Arry<T> &p) { return add(p, c); }


template <typename T>
Arry<T> substract(const Arry<T> &p1, const Arry<T> &p2)
{
    return add(p1, multiply(p2, T(-1.0)));
}

template <typename T>
Arry<T> substract(const Arry<T> &p, const T &c) { return add(p, -c); }

template <typename T>
Arry<T> substract(const T &c, const Arry<T> &p) { return substract(Arry<T>({c}), p); }


template <typename T>
Arry<T> multiply(const Arry<T> &p1, const Arry<T> &p2)
{
    CHECK_COEFS(p1, 1e-12);
    CHECK_COEFS(p2, 1e-12);

    // claculate the lenth of final polynomial
    const int &len = p1.size() + p2.size() - 1;

    // create an array holding result
    Arry<T> result(len, 0.0);

    for(unsigned int i=0; i<p1.size(); ++i)
    {
        const auto &c = p1[i];
        std::transform(p2.begin(), p2.end(), result.begin()+i, result.begin()+i,
            [&c](const T &x, const T &y)->T{return c*x+y;});
    }

    trim_leading_zeros(result, 1e-12);

    return result;
}

template <typename T>
Arry<T> multiply(const Arry<T> &p, const T &c)
{
    CHECK_COEFS(p, 1e-12);

    // special case. Note we use exactly zero here.
    if (std::abs(c) == 0.0) return {0.0};

    Arry<T> result(p);
    for(auto &it: result) it *= c;
    trim_leading_zeros(result, 1e-12);
    return result;
}

template <typename T>
Arry<T> multiply(const T &c, const Arry<T> &p) { return multiply(p, c); }


template <typename T>
Arry<T> divide(const Arry<T> &p1, const Arry<T> &p2, Arry<T> &r)
{
    CHECK_COEFS(p1, 1e-12);
    CHECK_COEFS(p2, 1e-12);

    // special case: p2 = constant
    if (p2.size() == 1)
    {
        r = Arry<T>(1, 0.0);
        return divide(p1, p2[0]);
    }

    // reset remainder's initial value to polynomial 1
    r = p1;

    // if polynomial 2 has a higher degree, then quotient is zero
    if (p1.size() < p2.size()) return Arry<T>(1, 0.0);

    // claculate the lenth of final polynomial
    const int &len = p1.size() - p2.size() + 1;

    // create vectors holding quotient
    Arry<T> Q(len);

    // alias of the leading coefficient of the polynomial 2
    const T &c = p2.back();

    for(int qi=len-1, ri=p1.size()-1; qi>=0; --qi, --ri)
    {
        Q[qi] = r[ri] / c;
        std::transform(p2.begin(), p2.end(), r.begin()+qi, r.begin()+qi,
            [&Q, &qi](T pp, T rr)->T{return rr-pp*Q[qi];});
    }

    r = Arry<T>(r.begin(), r.begin()+p2.size()-1);

    return Q;
}

template <typename T>
Arry<T> divide(const Arry<T> &p1, const Arry<T> &p2)
{
    Arry<T> R;
    return divide(p1, p2, R);
}

template <typename T>
Arry<T> divide(const Arry<T> &p, const T &c)
{
    CHECK_COEFS(p, 1e-12);

    // special case. Note we use exactly zero here.
    if (std::abs(c) == 0.0) throw exceptions::DivideByZero(__FL__);

    Arry<T> result(p);
    for(auto &it: result) it /= c;
    return result;
}


CArry to_CArry(const DArry &p)
{
    return CArry(p.begin(), p.end());
}

DArry to_DArry(const CArry &p, const double tol)
{
    DArry d(p.size());

# ifndef NDEBUG
    for(const auto &it: p) if (std::abs(it.imag()) > tol)
        throw exceptions::FoundComplexNumber(__FILE__, __LINE__, it);
# endif

    std::transform(p.begin(), p.end(), d.begin(),
            [](const Cmplx &x)->double{return x.real();});

    return d;
}


// GCD
template <typename T>
Arry<T> GCD(const Arry<T> &p1, const Arry<T> &p2, const double tol)
{
    CHECK_COEFS(p1, 1e-12);
    CHECK_COEFS(p2, 1e-12);

    Arry<T> a, b, q, r;
    double delta;

    a = p1;
    b = divide(p2, p2.back());
    auto f = [](const double &x, const T &y)->double{return x+std::norm(y);};

    int iter = 1;
    while(true)
    {
        q = divide(a, b, r);
        trim_leading_zeros(r, 1e-12);

        delta = std::sqrt(
                std::accumulate(std::begin(r), std::end(r), 0.0, f) /
                std::accumulate(std::begin(b), std::end(b), 0.0, f));

        if (delta < tol) break;

        a = b;
        b = divide(r, r.back());

        iter += 1;
        if (iter > 10000) throw exceptions::InfLoop(__FILE__, __LINE__);
    }

    return b;
}

// trim leading zero coefficients
template <typename T>
void trim_leading_zeros(Arry<T> &p, const double tol)
{
    while (p.size() > 1)
    {
        if (std::abs(p.back()) < tol) p.pop_back();
        else break;
    }
}


// find polynomial coefficients from roots
template <typename T>
Arry<T> to_coefficients(const T &l, const T* const &rts, const int len)
{
# ifndef NDEBUG
    if (len < 0) throw exceptions::NegativeCoeffsLength(__FL__, len);
# endif

    // polynomial f(x) = constant
    if (len == 0)
        return {l};
    else
        return multiply({-rts[0], 1.0}, to_coefficients(l, rts+1, len-1));
}

template <typename T>
Arry<T> to_coefficients(const T &l, const Arry<T> &rts)
{
    return to_coefficients(l, &rts[0], rts.size());
}


// explicit instantiation
template DArry add(const DArry &p1, const DArry &p2);
template CArry add(const CArry &p1, const CArry &p2);
template DArry add(const DArry &p, const double &c);
template CArry add(const CArry &p, const Cmplx &c);
template DArry add(const double &c, const DArry &p);
template CArry add(const Cmplx &c, const CArry &p);
template DArry substract(const DArry &p1, const DArry &p2);
template CArry substract(const CArry &p1, const CArry &p2);
template DArry substract(const DArry &p, const double &c);
template CArry substract(const CArry &p, const Cmplx &c);
template DArry substract(const double &c, const DArry &p);
template CArry substract(const Cmplx &c, const CArry &p);
template DArry multiply(const DArry &p1, const DArry &p2);
template CArry multiply(const CArry &p1, const CArry &p2);
template DArry multiply(const DArry &p, const double &c);
template CArry multiply(const CArry &p, const Cmplx &c);
template DArry multiply(const double &c, const DArry &p);
template CArry multiply(const Cmplx &c, const CArry &p);
template DArry divide(const DArry &p1, const DArry &p2, DArry &r);
template CArry divide(const CArry &p1, const CArry &p2, CArry &r);
template DArry divide(const DArry &p1, const DArry &p2);
template CArry divide(const CArry &p1, const CArry &p2);
template DArry divide(const DArry &p, const double &c);
template CArry divide(const CArry &p, const Cmplx  &c);

template DArry GCD(const DArry &p1, const DArry &p2, const double tol);
template CArry GCD(const CArry &p1, const CArry &p2, const double tol);

template void trim_leading_zeros(DArry &p, const double tol);
template void trim_leading_zeros(CArry &p, const double tol);

template DArry to_coefficients(const double &l, const double* const &rts, const int len);
template CArry to_coefficients(const Cmplx &l, const Cmplx* const &rts, const int len);
template DArry to_coefficients(const double &l, const DArry &rts);
template CArry to_coefficients(const Cmplx &l, const CArry &rts);


// operator << for DArry
std::ostream &operator<<(std::ostream &os, const simpoly::basic::DArry &v)
{
    if (v.size() == 0) return os;
    for(auto it=v.begin(); it<v.end()-1; ++it)
        os << *it << ", ";
    os << v.back();
    return os;
}

// operator << for CArry
std::ostream &operator<<(std::ostream &os, const simpoly::basic::CArry &v)
{
    if (v.size() == 0) return os;
    for(auto it=v.begin(); it<v.end()-1; ++it)
        os << *it << ", ";
    os << v.back();
    return os;
}

} // end of namespace basic
} // end of namespace simpoly
