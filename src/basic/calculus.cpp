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

    if (len == 1) return Arry<T>(1, 0.0);

    // initialize result using values in coeffs[1:]
    Arry<T> result(coeffs.begin()+1, coeffs.end());

    // alias to the beginning of `result` minus 1
    const auto &bg = result.data() - 1;

    // calculate correct values in-place
    std::for_each(result.begin(), result.end(),
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
    Arry<T> result(len+1, 0.0);

    // alias to the beginning of coeffs minus 1
    const auto &bg = coeffs.data() - 1;

    // perform result[i+1] = coeffs[i] / (i + 1)
    std::transform(coeffs.begin(), coeffs.end(), result.begin()+1,
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
