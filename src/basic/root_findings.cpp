/**
 * \file root_findings.cpp
 * \brief Root-finding functions of polynomials defined by std::valarray.
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
T newton_raphson(const Arry<T> &coeffs, const T guess, const double tol)
{
    T ans = guess;
    Arry<T> d = derivative(coeffs);
    
    for(unsigned iter=0; iter<10000; ++iter)
    {
        T value = evaluate(coeffs, ans),
          d_value = evaluate(d, ans),
          diff;
        
        if (value == 0.0) break; // found root; note we use exact 0.0 here
        
        if (d_value == 0.0) // to avoid devide by zero
        {
            ans *= 1.0001;
            value = evaluate(coeffs, ans);
            d_value = evaluate(d, ans);
        }
        diff = value / d_value;
        ans -= diff;
        if (std::abs(diff)/std::abs(((ans==0.0)?1.0:ans)) < tol) break;
        if (iter > 10000) throw exceptions::InfLoop(__FILE__, __LINE__);
    }
    
    return ans;
}


CArry aberth(const CArry &coeffs, const CArry &guess, const double tol)
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
    if (n == 0) return CArry(0);
    
    // initialize initial guess through copying
    CArry rts(guess);
    
    // create a vector indicating if each root has been found
    bool stop[n];
    std::fill(stop, stop+n, false);
    
    // derivative
    CArry d = derivative(coeffs);
    
    // an index to record the number of while iteration
    long iter = 0;
    
    while (std::any_of(stop, stop+n, [](bool b){return !b;}))
    {
        for(unsigned int i=0; i<n; ++i)
        {
            const Cmplx &zi = rts[i]; // alias
            Cmplx temp, delta;
            
            temp = evaluate(coeffs, zi) / evaluate(d, zi);
            
            delta = std::accumulate(
                std::begin(rts), std::begin(rts)+i, Cmplx(0.0, 0.0),
                [&zi] (const Cmplx &x, const Cmplx &y) -> Cmplx 
                    { return x + 1.0 / (zi - y); });
            
            delta += std::accumulate(
                std::begin(rts)+i+1, std::end(rts), Cmplx(0.0, 0.0),
                [&zi] (const Cmplx &x, const Cmplx &y) -> Cmplx 
                    { return x + 1.0 / (zi - y); });
            
            delta *= temp;
            delta = 1.0 - delta;
            delta = temp / delta;
            
            if (rts[i] == 0.0) stop[i] = true; // note we use exact zero here
            if ((std::abs(delta)/std::abs(rts[i])) < tol) stop[i] = true;
            
            rts[i] -= delta;
        }
        
        // check the number of iterations
        iter += 1;
        if (iter > 10000) throw exceptions::InfLoop(__FILE__, __LINE__);
    }
    
    return rts;
}


CArry aberth(const CArry &coeffs, const DArry &guess, const double tol)
{
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return CArry(0);
    
    // make a copy of guess with complex type, and perturbation in imag
    CArry G = to_CArry(guess);
    
    // add perturbation to imaginary parts
    std::for_each(std::begin(G), std::end(G), [&tol](Cmplx &i){i.imag(tol);});
    
    return aberth(coeffs, G, tol);
}


CArry aberth(const DArry &coeffs, const CArry &guess, const double tol)
{
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return CArry(0);
    
    return aberth(to_CArry(coeffs), guess, tol);
}


CArry aberth(const DArry &coeffs, const DArry &guess, const double tol)
{
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return CArry(0);
    
    // make a copy of guess with complex type
    CArry G = to_CArry(guess);
    
    // add perturbation to imaginary parts
    std::for_each(std::begin(G), std::end(G), [&tol](Cmplx &i){i.imag(tol);});
    
    return aberth(coeffs, G, tol);
}


CArry aberth(const CArry &coeffs, const double tol)
{
# ifndef NDEBUG
    using namespace exceptions;
    if (coeffs.size() == 0) throw ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return CArry(0);
    
    // initialize initial guess through copying
    CArry guess(coeffs.size()-1);
    
    std::iota(std::begin(guess), std::end(guess), 0.0);
    guess = std::pow(Cmplx(0.5, 0.5), guess);
    
    return aberth(coeffs, guess, tol);
}


CArry aberth(const DArry &coeffs, const double tol)
{
# ifndef NDEBUG
    using namespace exceptions;
    if (coeffs.size() == 0) throw ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return CArry(0);
    
    // initialize initial guess through copying
    CArry guess(coeffs.size()-1);
    
    std::iota(std::begin(guess), std::end(guess), 0.0);
    guess = std::pow(Cmplx(0.5, 0.5), guess);
    
    return aberth(coeffs, guess, tol);
}


DArry aberth_real(const DArry &coeffs, const DArry &guess, 
        const double tol, const bool no_ignore_cmplx)
{
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return DArry(0);
    
    // get roots using algorithms that is specifically for complex roots
    CArry z = aberth(coeffs, guess, tol);;
    
    // create an array to hold real roots
    DArry roots(z.size());
    
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


DArry aberth_real(const DArry &coeffs, 
        const double tol, const bool no_ignore_cmplx)
{
# ifndef NDEBUG
    using namespace exceptions;
    if (coeffs.size() == 0) throw ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // if degree is 0, no root exists; return a zero-length array
    if (coeffs.size() == 1) return DArry(0);
    
    // initialize initial guess through copying
    CArry guess(coeffs.size()-1);
    
    std::iota(std::begin(guess), std::end(guess), 0.0);
    guess = std::pow(Cmplx(0.5, 0.5), guess);
    
    // get roots
    CArry z = aberth(coeffs, guess, tol);;
    
    // create an array to hold real roots
    DArry roots(z.size());
    
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


// explicit instantiation
template double newton_raphson(const DArry &coeffs, const double guess, const double tol);
template Cmplx newton_raphson(const CArry &coeffs, const Cmplx guess, const double tol);

} // end of namespace basic
} // end of namespace simpoly
