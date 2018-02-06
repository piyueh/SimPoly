/**
 * \file calculus.cpp
 * \brief Calculus of polynomials defined by a std::valarray.
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
Arry<T> derivative(const Arry<T> &coeffs)
{
    // alias to the length of provided coefficient array
    const auto &len = coeffs.size();
    
# ifndef NDEBUG
    if (len == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    if (len == 1) return Arry<T>(0.0, 1);
    
    // initialize result using values in coeffs[1:]
    Arry<T> result(coeffs[std::slice(1, len-1, 1)]);
    
    // alias to the beginning of `result` minus 1
    const auto &bg = std::begin(result) - 1;
    
    // calculate correct values in-place
    std::for_each(std::begin(result), std::end(result),
            [&bg](T &x){x *= T(&x - bg);});
    
    return result;
}


template <typename T>
Arry<T> integral(const Arry<T> &coeffs)
{
    // alias to the length of provided coefficient array
    const auto &len = coeffs.size();
    
# ifndef NDEBUG
    if (len == 0) throw exceptions::ZeroCoeffsLength(__FILE__, __LINE__);
# endif
    
    // initialize result with len+1 zeros
    Arry<T> result(0.0, len+1);
    
    // alias to the beginning of coeffs minus 1
    const auto &bg = std::begin(coeffs) - 1;
    
    // perform result[i+1] = coeffs[i] / (i + 1)
    std::transform(std::begin(coeffs), std::end(coeffs), std::begin(result)+1,
            [&bg](const T &x)->T{return x / T(&x - bg);});
    
    return result;
}


// explicit instantiation
template DArry derivative(const DArry &coeffs);
template CArry derivative(const CArry &coeffs);
template DArry integral(const DArry &coeffs);
template CArry integral(const CArry &coeffs);
} // end of namespace basic
} // end of namespace simpoly
