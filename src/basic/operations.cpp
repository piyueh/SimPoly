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
    Arry<T> result;
    if (p1.size() >= p2.size())
    {
        result.resize(p1.size());
        std::copy(std::begin(p2), std::end(p2), std::begin(result));
        result += p1;
    }
    else
    {
        result.resize(p2.size());
        std::copy(std::begin(p1), std::end(p1), std::begin(result));
        result += p2;
    }
    
    return result;
}
    

template <typename T>
Arry<T> multiply(const Arry<T> &p1, const Arry<T> &p2)
{
# ifndef NDEBUG
    if (p1.size() == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
    if (p2.size() == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // claculate the lenth of final polynomial
    const int &len = p1.size() + p2.size() - 1;
    
    // create an array holding result
    Arry<T> result(0.0, len);
    
    for(unsigned int i=0; i<p1.size(); ++i)
    {
        const auto &c = p1[i];
        std::transform(
            std::begin(p2), std::end(p2), 
            std::begin(result)+i, std::begin(result)+i,
            [&c](const T &x, const T &y)->T{return c*x+y;});
    }
    
    return result;
}
    

template <typename T>
Arry<T> divide(const Arry<T> &p1, const Arry<T> &p2, Arry<T> &r)
{
# ifndef NDEBUG
    if (p1.size() == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
    if (p2.size() == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // reset remainder's initial value to polynomial 1
    r = p1;
    
    // if polynomial 2 has a higher degree, then quotient is zero
    if (p1.size() < p2.size()) return Arry<T>({0.0});
    
    // claculate the lenth of final polynomial
    const int &len = p1.size() - p2.size() + 1;
    
    // create vectors holding quotient
    Arry<T> Q(len);
    
    // alias of the leading coefficient of the polynomial 2
    const T &c = *(std::end(p2)-1);
    
    for(int qi=len-1, ri=p1.size()-1; qi>=0; --qi, --ri)
    {
        Q[qi] = r[ri] / c;
        std::transform(
            std::begin(p2), std::end(p2), std::begin(r)+qi, std::begin(r)+qi,
            [&Q, &qi](T pp, T rr)->T{return rr-pp*Q[qi];});
    }
    
    r = Arry<T>(r[std::slice(0, p2.size()-1, 1)]);
    
    return Q;
}


template <typename T>
Arry<T> divide(const Arry<T> &p1, const Arry<T> &p2)
{
    Arry<T> R;
    return divide(p1, p2, R);
}


CArry to_CArry(const DArry &p)
{
    CArry result(p.size());
    std::copy(std::begin(p), std::end(p), std::begin(result));
    return result;
}


// explicit instantiation
template DArry add(const DArry &p1, const DArry &p2);
template CArry add(const CArry &p1, const CArry &p2);
template DArry multiply(const DArry &p1, const DArry &p2);
template CArry multiply(const CArry &p1, const CArry &p2);
template DArry divide(const DArry &p1, const DArry &p2, DArry &r);
template CArry divide(const CArry &p1, const CArry &p2, CArry &r);
template DArry divide(const DArry &p1, const DArry &p2);
template CArry divide(const CArry &p1, const CArry &p2);


// GCD
template <typename T>
Arry<T> GCD(const Arry<T> &p1, const Arry<T> &p2, const double tol)
{
    Arry<T> a, b, q, r;
    double delta;
    
    a = p1;
    b = p2 / (*(std::end(p2)-1));
    auto f = [](const double &x, const T &y) -> double { return x + std::norm(y); };
    
    int iter = 1;
    do
    {
        q = divide(a, b, r);
        
        delta = std::sqrt(
                std::accumulate(std::begin(r), std::end(r), 0.0, f) /
                std::accumulate(std::begin(b), std::end(b), 0.0, f));
        
        a = b;
        b = r / (*(std::end(r)-1));
        
        iter += 1;
        if (iter > 10000) throw exceptions::InfLoop(__FILE__, __LINE__);
    } while(delta > tol);
    
    return a;
}

template DArry GCD(const DArry &p1, const DArry &p2, const double tol);
template CArry GCD(const CArry &p1, const CArry &p2, const double tol);

} // end of namespace basic
} // end of namespace simpoly


// operator << for DArry
std::ostream &operator<<(std::ostream &os, const simpoly::basic::DArry &v)
{
    for(auto it=std::begin(v); it<std::end(v)-1; ++it)
        os << *it << ", ";
    os << (*(std::end(v)-1));
    return os;
}


// operator << for CArry
std::ostream &operator<<(std::ostream &os, const simpoly::basic::CArry &v)
{
    for(auto it=std::begin(v); it<std::end(v)-1; ++it)
        os << *it << ", ";
    os << (*(std::end(v)-1));
    return os;
}
