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
    
# ifndef NDEBUG
    using namespace exceptions;
    if (len == 0) throw ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // if degree is 0, no root exists; return a zero-length array
    if (len == 1) return CArry(0);
    
    // initialize initial guess through copying
    CArry rts(guess);
    
    // create a vector indicating if each root has been found
    bool stop[guess.size()];
    std::fill(stop, stop+guess.size(), false);
    
    // derivative
    CArry d = derivative(coeffs);
    
    // an index to record the number of while iteration
    long iter = 0;
    
    while (std::any_of(stop, stop+guess.size(), [](bool b){return !b;}))
    {
        for(unsigned int i=0; i<guess.size(); ++i)
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


CArry yan_and_chieng_2006(const CArry &coeffs, const double tol)
{
# ifndef NDEBUG
    if (coeffs.size() == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
# endif

    CArry     drv = derivative(coeffs); // direvative
    CArry     agcd = GCD(coeffs, drv); // approximated GCD
    CArry     q_coeffs = divide(coeffs, agcd); // f(x) / GCD
    CArry     q_drv = divide(drv, agcd); // direvative / GCD
    CArry     drv_q_coeffs = derivative(q_coeffs);
    
    // simple roots from q_coeffs
    CArry       simples1 = aberth(q_coeffs, tol); 
    
    // refine simple roots with original polynomial `coeffs`
    CArry       simples2 = aberth(coeffs, simples1, tol);
    
    // an array for final results
    CArry       result(coeffs.size() - 1);
    
    unsigned k = 0; // index for result
    for(unsigned i=0; i<simples1.size(); ++i)
    {
        // claculate multiplicity
        unsigned m = int((evaluate(q_drv, simples1[i]) /
                evaluate(drv_q_coeffs, simples1[i])).real() + 0.5);
        
        // choose the best one
        if (m > 1)
        {
            double e1 = std::norm(evaluate(drv, simples1[i])) + 
                        std::norm(evaluate(coeffs, simples1[i])),
                   e2 = std::abs(evaluate(drv, simples2[i])) + 
                        std::norm(evaluate(coeffs, simples2[i]));
            
            if (e1 < e2) result[k] = simples1[i];
            else result[k] = simples2[i];
            k += 1;
        
            // duplicate multiple roots
            for(unsigned mi=1; mi<m; ++mi)
            {
                result[k] = result[k-1];
                k += 1;
            }
        }
        else
        {
            double e1 = std::abs(evaluate(coeffs, simples1[i])),
                   e2 = std::abs(evaluate(coeffs, simples2[i]));
            
            if (e1 < e2) result[k] = simples1[i];
            else result[k] = simples2[i];
            k += 1;
        }
    }
    
    return result;
}


CArry yan_and_chieng_2006(const DArry &coeffs, const double tol)
{
    CArry C = to_CArry(coeffs);
    return yan_and_chieng_2006(C, tol);
}


// explicit instantiation
template double newton_raphson(const DArry &coeffs, const double guess, const double tol);
template Cmplx newton_raphson(const CArry &coeffs, const Cmplx guess, const double tol);

} // end of namespace basic
} // end of namespace simpoly
