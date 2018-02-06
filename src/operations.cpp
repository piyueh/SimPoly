/**
 * \file operations.cpp
 * \brief Operations of polynomials defined by a std::vector<double>.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-29
 */


# include <iostream>
# include <complex>
# include <vector>
# include <valarray>
# include <algorithm>
# include <numeric>
# include <cmath>

# include "exceptions.h"


namespace simpoly
{
namespace op
{
    

std::valarray<double> add(const std::valarray<double> &p1, 
        const std::valarray<double> &p2)
{
    std::valarray<double> result;
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
    

std::valarray<double> multiply(const std::valarray<double> &p1, 
        const std::valarray<double> &p2)
{
# ifndef NDEBUG
    if (p1.size() == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
    if (p2.size() == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // claculate the lenth of final polynomial
    const int &len = p1.size() + p2.size() - 1;
    
    // create an array holding result
    std::valarray<double> result(0.0, len);
    
    for(unsigned int i=0; i<p1.size(); ++i)
    {
        const auto &c = p1[i];
        std::transform(std::begin(p2), std::end(p2), 
                std::begin(result)+i, std::begin(result)+i,
                [&c](const double &x, const double &y)->double{return c*x+y;});
    }
    
    return result;
}
    

std::valarray<double> divide(const std::valarray<double> &p1, 
        const std::valarray<double> &p2, std::valarray<double> &r)
{
# ifndef NDEBUG
    if (p1.size() == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
    if (p2.size() == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // reset remainder's initial value to polynomial 1
    r = p1;
    
    // if polynomial 2 has a higher degree, then quotient is zero
    if (p1.size() < p2.size()) return std::valarray<double>({0.0});
    
    // claculate the lenth of final polynomial
    const int &len = p1.size() - p2.size() + 1;
    
    // create vectors holding quotient
    std::valarray<double> Q(len);
    
    // alias
    const double &c = *(std::end(p2)-1);
    
    int i = Q.size();
    for(int qi=len-1, ri=p1.size()-1; qi>=0; --qi, --ri)
    {
        Q[qi] = r[ri] / c;
        std::transform(
                std::begin(p2), std::end(p2), std::begin(r)+qi, std::begin(r)+qi,
                [&Q, &qi](double pp, double rr)->double{return rr-pp*Q[qi];});
    }
    
    r = std::valarray<double>(r[std::slice(0, p2.size()-1, 1)]);
    
    return Q;
}


std::valarray<double> divide(
        const std::valarray<double> &p1, const std::valarray<double> &p2)
{
    std::valarray<double> R;
    return divide(p1, p2, R);
}

    
double evaluate(const double* const &bg, const int len, const double x)
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


double evaluate(const std::valarray<double> &coeffs, const double x)
{
    return evaluate(&coeffs[0], coeffs.size(), x);
}
    
    
std::complex<double> evaluate(const std::complex<double>* const &bg, 
        const int len, const std::complex<double> x)
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


std::complex<double> evaluate(const std::valarray<std::complex<double>> &coeffs, 
        const std::complex<double> x)
{
    return evaluate(&coeffs[0], coeffs.size(), x);
}


double evaluate_root(const double l, 
            const double* const &bg, const int degree, const double x)
{
# ifndef NDEBUG
    if (degree < 0) throw exceptions::NegativeDegree(__FILE__, __LINE__, degree);
# endif
    
    if (degree == 0)
        return l;
    else
        return (x - *bg) * evaluate_root(l, bg+1, degree-1, x);
}


double evaluate_root(const double l, 
            const std::valarray<double> &roots, const double x)
{
    return evaluate_root(l, &roots[0], roots.size(), x);
}


std::complex<double> evaluate_root(
        const double l, const std::complex<double>* const &bg, 
        const int degree, const std::complex<double> x)
{
# ifndef NDEBUG
    if (degree < 0) throw exceptions::NegativeDegree(__FILE__, __LINE__, degree);
# endif
    
    if (degree == 0)
        return l;
    else
        return (x - *bg) * evaluate_root(l, bg+1, degree-1, x);
}


std::complex<double> evaluate_root(const double l, 
        const std::valarray<std::complex<double>> &roots, 
        const std::complex<double> x)
{
    return evaluate_root(l, &roots[0], roots.size(), x);
}


std::valarray<double> derivative(const std::valarray<double> &coeffs)
{
    // alias to the length of provided coefficient array
    const auto &len = coeffs.size();
    
# ifndef NDEBUG
    if (len == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    if (len == 1)
        return std::valarray<double>(0.0, 1);
    
    // initialize result using values in coeffs[1:]
    std::valarray<double> result(coeffs[std::slice(1, len-1, 1)]);
    
    // alias to the beginning of result minus 1
    const auto &bg = std::begin(result) - 1;
    
    // calculate correct values in-place
    std::for_each(std::begin(result), std::end(result),
            [&bg](double &x){x *= (&x - bg);});
    
    return result;
}


std::valarray<std::complex<double>> derivative(
        const std::valarray<std::complex<double>> &coeffs)
{
    // alias to the length of provided coefficient array
    const auto &len = coeffs.size();
    
# ifndef NDEBUG
    if (len == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    if (len == 1)
        return std::valarray<std::complex<double>>(0.0, 1);
    
    // initialize result using values in coeffs[1:]
    std::valarray<std::complex<double>> result(coeffs[std::slice(1, len-1, 1)]);
    
    // alias to the beginning of result minus 1
    const auto &bg = std::begin(result) - 1;
    
    // calculate correct values in-place
    std::for_each(std::begin(result), std::end(result),
            [&bg](std::complex<double> &x){x *= (&x - bg);});
    
    return result;
}


std::valarray<double> integral(const std::valarray<double> &coeffs)
{
    // alias to the length of provided coefficient array
    const auto &len = coeffs.size();
    
# ifndef NDEBUG
    if (len == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // initialize result with len+1 zeros
    std::valarray<double> result(0.0, len+1);
    
    // alias to the beginning of coeffs minus 1
    const auto &bg = std::begin(coeffs) - 1;
    
    // perform result[i+1] = coeffs[i] / (i + 1)
    std::transform(std::begin(coeffs), std::end(coeffs), std::begin(result)+1,
            [&bg](const double &x)->double{return x / (&x - bg);});
    
    return result;
}


double find_root_NRM(const std::valarray<double> &coeffs,
        const double guess, const double tol)
{
    double ans, diff;
    std::valarray<double> d = derivative(coeffs);
    
    ans = guess;
    
    for(unsigned iter=0; iter<10000; ++iter)
    {
        double value = evaluate(coeffs, ans),
               d_value = evaluate(d, ans);
        
        if (d_value == 0.0)
        {
            if (value == 0.0) break; // found multiple root
            else
            {
                ans *= 1.01;
                value = evaluate(coeffs, ans);
                d_value = evaluate(d, ans);
            }
        }
        diff = value / d_value;
        ans -= diff;
        if (std::abs(diff/((ans==0.0)?1.0:ans)) < tol) break;
        if (iter > 10000) throw exceptions::InfLoop(__FILE__, __LINE__);
    }
    
    return ans;
}


std::complex<double> find_root_NRM(
        const std::valarray<std::complex<double>> &coeffs,
        const std::complex<double> guess, const double tol)
{
    std::complex<double> ans, diff;
    std::valarray<std::complex<double>> d = derivative(coeffs);
    
    ans = guess;
    
    for(unsigned iter=0; iter<10000; ++iter)
    {
        std::complex<double> 
            value = evaluate(coeffs, ans),
            d_value = evaluate(d, ans);
        
        if (d_value == 0.0)
        {
            if (value == 0.0) break; // found multiple root
            else
            {
                ans *= 1.01;
                value = evaluate(coeffs, ans);
                d_value = evaluate(d, ans);
            }
        }
        diff = value / d_value;
        ans -= diff;
        if ((std::abs(diff)/std::abs((ans==0.0)?1.0:ans)) < tol) break;
        if (iter > 10000) throw exceptions::InfLoop(__FILE__, __LINE__);
    }
    
    return ans;
}


std::valarray<std::complex<double>> find_roots_complex(
        const std::valarray<std::complex<double>> &coeffs, 
        const std::valarray<std::complex<double>> &guess, const double tol)
{
    // alias to the length of provided coefficient array
    const auto &len = coeffs.size();
    
    // expected number of roots
    const int &n = len - 1;
    
# ifndef NDEBUG
    using namespace exceptions;
    if (len == 0) throw ZeroCoeffsLength(__FILE__, __LINE__);
    if (n != guess.size()) throw UnmatchedLength(__FILE__, __LINE__, n, guess.size());
# endif
    
    // if degree is 0, no root exists; return a zero-length array
    if (n == 0) return std::valarray<std::complex<double>>(0);
    
    // initialize initial guess through copying
    std::valarray<std::complex<double>> roots(guess);
    
    // create a vector indicating if each root has been found
    std::vector<bool> stop(n, false);
    
    // derivative
    auto d = derivative(coeffs);
    
    // an index to record the number of where iteration
    long iter = 0;
    
    while (std::any_of(stop.begin(), stop.end(), [](bool b){return !b;}))
    {
        for(unsigned int i=0; i<n; ++i)
        {
            const std::complex<double> &zi = roots[i];
            std::complex<double> temp, delta;
            
            temp = evaluate(coeffs, zi) / evaluate(d, zi);
            
            delta = std::accumulate(
                std::begin(roots), std::begin(roots)+i, std::complex<double>(0.0, 0.0),
                [&zi] (const std::complex<double> &x, const std::complex<double> &y)
                -> std::complex<double> { return x + 1.0 / (zi - y); });
            
            delta += std::accumulate(
                std::begin(roots)+i+1, std::end(roots), std::complex<double>(0.0, 0.0),
                [&zi] (const std::complex<double> &x, const std::complex<double> &y)
                -> std::complex<double> { return x + 1.0 / (zi - y); });
            
            delta *= temp;
            delta = 1.0 - delta;
            delta = temp / delta;
            
            if ((roots[i].real() == 0.0) && 
                    (roots[i].imag() == 0.0)) stop[i] = true;
            if ((std::abs(delta)/std::abs(roots[i])) < tol) stop[i] = true;
            
            roots[i] -= delta;
        }
        
        // check the number of iterations
        iter += 1;
        if (iter > 10000) throw exceptions::InfLoop(__FILE__, __LINE__);
    }
    
    return roots;
}


std::valarray<std::complex<double>> find_roots_complex(
        const std::valarray<std::complex<double>> &coeffs, 
        const std::valarray<double> &guess, const double tol)
{
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return std::valarray<std::complex<double>>(0);
    
    // make a copy of guess with complex type, and perturbation in imag
    std::valarray<std::complex<double>> G(guess.size());
    std::copy(std::begin(guess), std::end(guess), std::begin(G));
    std::for_each(std::begin(G), std::end(G),
            [](std::complex<double> &i){i.imag(1e-8);});
    
    return find_roots_complex(coeffs, G, tol);
}


std::valarray<std::complex<double>> find_roots_complex(
        const std::valarray<double> &coeffs, 
        const std::valarray<std::complex<double>> &guess, const double tol)
{
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return std::valarray<std::complex<double>>(0);
    
    // make a copy of coeffs with complex type
    std::valarray<std::complex<double>> C(coeffs.size());
    std::copy(std::begin(coeffs), std::end(coeffs), std::begin(C));
    
    return find_roots_complex(C, guess, tol);
}


std::valarray<std::complex<double>> find_roots_complex(
        const std::valarray<double> &coeffs, 
        const std::valarray<double> &guess, const double tol)
{
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return std::valarray<std::complex<double>>(0);
    
    // make a copy of guess with complex type
    std::valarray<std::complex<double>> G(guess.size());
    std::copy(std::begin(guess), std::end(guess), std::begin(G));
    std::for_each(std::begin(G), std::end(G),
            [](std::complex<double> &i){i.imag(1e-8);});
    
    return find_roots_complex(coeffs, G, tol);
}


std::valarray<std::complex<double>> find_roots_complex(
        const std::valarray<std::complex<double>> &coeffs, const double tol)
{
# ifndef NDEBUG
    using namespace exceptions;
    if (coeffs.size() == 0) throw ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return std::valarray<std::complex<double>>(0);
    
    // initialize initial guess through copying
    std::valarray<std::complex<double>> guess(coeffs.size()-1);
    
    std::iota(std::begin(guess), std::end(guess), 0.0);
    guess = std::pow(std::complex<double>(0.5, 0.5), guess);
    
    return find_roots_complex(coeffs, guess, tol);
}


std::valarray<std::complex<double>> find_roots_complex(
        const std::valarray<double> &coeffs, const double tol)
{
# ifndef NDEBUG
    using namespace exceptions;
    if (coeffs.size() == 0) throw ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return std::valarray<std::complex<double>>(0);
    
    // initialize initial guess through copying
    std::valarray<std::complex<double>> guess(coeffs.size()-1);
    
    std::iota(std::begin(guess), std::end(guess), 0.0);
    guess = std::pow(std::complex<double>(0.5, 0.5), guess);
    
    return find_roots_complex(coeffs, guess, tol);
}


std::valarray<double> find_roots(
        const std::valarray<double> &coeffs, const std::valarray<double> &guess, 
        const double tol, const bool no_ignore_cmplx)
{
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return std::valarray<double>(0);
    
    // get roots using algorithms that is specifically for complex roots
    std::valarray<std::complex<double>> z = find_roots_complex(coeffs, guess, tol);;
    
    // create an array to hold real roots
    std::valarray<double> roots(z.size());
    
    // eliminate imagine part
    for(unsigned int i=0; i<z.size(); ++i)
    {
        if ((no_ignore_cmplx) && (std::abs(z[i].imag()) >= tol))
            throw exceptions::ComplexRoot(
                    __FILE__, __LINE__, z[i].real(), z[i].imag());
        
        roots[i] = z[i].real();
    }
    
    return roots;
}


std::valarray<double> find_roots(
        const std::valarray<double> &coeffs, 
        const double tol, const bool no_ignore_cmplx)
{
# ifndef NDEBUG
    using namespace exceptions;
    if (coeffs.size() == 0) throw ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return std::valarray<double>(0);
    
    // initialize initial guess through copying
    std::valarray<std::complex<double>> guess(coeffs.size()-1);
    
    std::iota(std::begin(guess), std::end(guess), 0.0);
    guess = std::pow(std::complex<double>(0.5, 0.5), guess);
    
    // get roots using algorithms that is specifically for complex roots
    std::valarray<std::complex<double>> z = find_roots_complex(coeffs, guess, tol);;
    
    // create an array to hold real roots
    std::valarray<double> roots(z.size());
    
    // eliminate imagine part
    for(unsigned int i=0; i<z.size(); ++i)
    {
        if ((no_ignore_cmplx) && (std::abs(z[i].imag()) >= tol))
            throw exceptions::ComplexRoot(
                    __FILE__, __LINE__, z[i].real(), z[i].imag());
        
        roots[i] = z[i].real();
    }
    
    return roots;
}

} // end of namespace op
} // end of namespace simpoly
