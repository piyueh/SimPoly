/**
 * \file basic.h
 * \brief Basic operations of polynomials defined by a std::valarray.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-02-06
 */

# pragma once

# include <complex>
# include <valarray>


namespace simpoly
{
namespace basic
{

/** \brief Alias of std::valarray<...>.
 *
 * \tparam T Basic type of each entry in the std::valarray. */
template <typename T> using Arry = std::valarray<T>;

/** \brief Alias of std::complex<double>. */
typedef std::complex<double> Cmplx;
    
/** \brief Alias of std::valarray<double>. */
typedef Arry<double> DArry;

/** \brief Alias of std::valarray<std::complex<double>>. */
typedef Arry<Cmplx> CArry;

/**
 * \brief Convert a DArry to CArry.
 *
 * \param p [in] A DArry (i.e., std::valarray<double>).
 *
 * \return A CArry (i.e., std::valarray<std::complex<double>>).
 */
CArry to_CArry(const DArry &p);
    
    
/**
 * \brief Addition of two polynomials defined by std::valarray.
 *
 * \tparam T Basic type of each entry in the std::valarray.
 * \param p1 [in] Coefficients of the first polynomial.
 * \param p2 [in] Coefficients of the second polynomial.
 *
 * \return Resulting polynomial.
 */
template <typename T>
Arry<T> add(const Arry<T> &p1, const Arry<T> &p2);

/**
 * \brief Multiplication of two polynomials defined by std::valarray.
 *
 * \tparam T Basic type of each entry in the std::valarray.
 * \param p1 [in] Coefficients of the first polynomial.
 * \param p2 [in] Coefficients of the second polynomial.
 *
 * \return Resulting polynomial.
 */
template <typename T>
Arry<T> multiply(const Arry<T> &p1, const Arry<T> &p2);

/**
 * \brief Division of two polynomials defined by std::valarray.
 *
 * \tparam T Basic type of each entry in the std::valarray.
 * \param p1 [in] Coefficients of the first polynomial.
 * \param p2 [in] Coefficients of the second polynomial.
 * \param r [out] Remainder polynomial.
 *
 * \return Quotient polynomial.
 * 
 * This function calculates the Q(x) and r(x) in p1(x) = p2(x) * Q(x) + r(x).
 */
template <typename T>
Arry<T> divide(const Arry<T> &p1, const Arry<T> &p2, Arry<T> &r);

/**
 * \brief Division of two polynomials defined by std::valarray.
 *
 * \tparam T Basic type of each entry in the std::valarray.
 * \param p1 [in] Coefficients of the first polynomial.
 * \param p2 [in] Coefficients of the second polynomial.
 *
 * \return Quotient polynomial.
 * 
 * This function calculates the Q(x) in p1(x) = p2(x) * Q(x). And this function
 * assumes that there is no remainder polynomial, i.e., r(x) = 0. If there is
 * actually a non-trivial remainder, this function will just discard the 
 * remainder and won't check or return any exception. So it's users'
 * responsibility to be sure there won't be any remainder.
 */
template <typename T>
Arry<T> divide(const Arry<T> &p1, const Arry<T> &p2);

/**
 * \brief Find greatest common divisor of two polynomials.
 *
 * \tparam T Base type of the entries in coefficient array.
 * \param p1 [in] A std::valarray representing the first polynomial.
 * \param p2 [in] A std::valarray representing the second polynomial.
 * \param tol [in] Tolerance (default 1e-8).
 *
 * \return The greatest common divisor polynomial.
 */
template <typename T>
Arry<T> GCD(const Arry<T> &p1, const Arry<T> &p2, const double tol=1e-8);


/**
 * \brief Evaluate polynomial value at a specified location.
 * 
 * This function evaluates the polynomial value at location x. The coefficient
 * array of the polynomial is defined by a raw pointer.
 *
 * \tparam T Basic type of each entry in the coefficient array.
 * \param bg [in] Pointer to the first entry in coefficient array.
 * \param len [in] The length of coefficient array.
 * \param x [in] The specified location.
 *
 * \return Polynomial value.
 */
template <typename T> 
T evaluate(const T* const &bg, const int len, const T x);

/**
 * \brief Evaluate polynomial value at a specified location.
 * 
 * This function evaluates the polynomial value at location x. The coefficient
 * array of the polynomial is defined by std::valarray.
 *
 * \tparam T Basic type of each entry in the std::valarray.
 * \param bg [in] Pointer to the first entry in coefficient array.
 * \param len [in] The length of coefficient array.
 * \param x [in] The specified location.
 *
 * \return Polynomial value.
 */
template <typename T>
T evaluate(const Arry<T> &coeffs, const T x);

/**
 * \brief Evaluate polynomial value at a specified location by providing roots.
 * 
 * This function evaluates the polynomial value at location x. The polynomial is
 * defined by a raw array of roots and a leading coefficient.
 *
 * \tparam T Basic type of each entry in the root array.
 * \param l [in] The coefficient of the highest-degree term.
 * \param bg [in] Pointer to the first entry in root array.
 * \param degree [in] The length of root array.
 * \param x [in] The specified location.
 *
 * \return Polynomial value.
 */
template <typename T>
T evaluate_from_root(const T l, const T* const &bg, const int degree, const T x);

/**
 * \brief Evaluate polynomial value at a specified location by providing roots.
 * 
 * This function evaluates the polynomial value at location x. The polynomial is
 * defined by a std::valarray of roots and a leading coefficient.
 *
 * \tparam T Basic type of each entry in the root std::valarray.
 * \param l [in] The coefficient of the highest-degree term.
 * \param roots [in] A std::valarray of roots.
 * \param x [in] The specified location.
 *
 * \return Polynomial value.
 */
template <typename T>
T evaluate_from_root(const T l, const Arry<T> &roots, const T x);


/**
 * \brief Obtain the derived polynomial of a polynomial.
 *
 * \tparam T Basic type of each entry in the root std::valarray.
 * \param coeffs [in] A polynomial defined by std::valarray.
 *
 * \return Derived polynomial.
 */
template <typename T>
Arry<T> derivative(const Arry<T> &coeffs);

/**
 * \brief Obtain the undeterministic integral of a polynomial.
 *
 * \tparam T Basic type of each entry in the root std::valarray.
 * \param coeffs [in] A polynomial defined by std::valarray.
 *
 * \return Undeterministic integral.
 */
template <typename T>
Arry<T> integral(const Arry<T> &coeffs);

    
/**
 * \brief Basic Newton-Raphson root-finding method.
 * 
 * Newton-Raphson method can only find one root each time and require a close
 * initial guess to the target root. This is the basic version, not the modified
 * Newton-Raphson, so the convergence of a multiple root is difficult.
 *
 * \tparam T The basic data type of the coefficient array.
 * \param coeffs [in] A std::valarray representing polynomial coefficients.
 * \param guess [in] An initial guess.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return Root.
 */
template <typename T>
T newton_raphson(const Arry<T> &coeffs, const T guess, const double tol=1e-13);

/**
 * \brief Aberth method for root finding.
 * 
 * Aberth method can find all roots at once. However, it only works on complex
 * plane. So the returned roots are complex numbers, even the imaginary part is
 * very close to zero (or even equal to zero).
 * 
 * Aberth method is good for simple roots but bad when there are multiple roots.
 *
 * \param coeffs [in] A std::valarray representing polynomial coefficients.
 * \param guess [in] A std::valarray of initial guess to all roots.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::valarray of all roots.
 */
CArry aberth(const CArry &coeffs, const CArry &guess, const double tol=1e-13);

/**
 * \brief Aberth method for root finding.
 * 
 * An overloaded version that provided initial guess is a std::valarray of
 * real numbers. The output is still a std::valarray of complex numbers.
 *
 * \param coeffs [in] A std::valarray representing polynomial coefficients.
 * \param guess [in] A std::valarray of initial guess to all roots.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::valarray of all roots.
 */
CArry aberth(const CArry &coeffs, const DArry &guess, const double tol=1e-13);

/**
 * \brief Aberth method for root finding.
 * 
 * An overloaded version that provided coefficients is a std::valarray of
 * real numbers. The output is still a std::valarray of complex numbers.
 *
 * \param coeffs [in] A std::valarray representing polynomial coefficients.
 * \param guess [in] A std::valarray of initial guess to all roots.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::valarray of all roots.
 */
CArry aberth(const DArry &coeffs, const CArry &guess, const double tol=1e-13);

/**
 * \brief Aberth method for root finding.
 * 
 * An overloaded version that both provided coefficients and initial guess are 
 * `std::valarray`s of real numbers. The output is still a std::valarray of 
 * complex numbers.
 *
 * \param coeffs [in] A std::valarray representing polynomial coefficients.
 * \param guess [in] A std::valarray of initial guess to all roots.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::valarray of all roots.
 */
CArry aberth(const DArry &coeffs, const DArry &guess, const double tol=1e-13);

/**
 * \brief Aberth method for root finding.
 * 
 * An overloaded version of Aberth method that uses default initial guess,
 * instead of user-provided initial guess.
 *
 * \param coeffs [in] A std::valarray representing polynomial coefficients.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::valarray of all roots.
 */
CArry aberth(const CArry &coeffs, const double tol=1e-13);

/**
 * \brief Aberth method for root finding.
 * 
 * An overloaded version of Aberth method that uses default initial guess,
 * instead of user-provided initial guess. Also, the coefficient array is a
 * std::valarray of real numbers, not complex numbers.
 *
 * \param coeffs [in] A std::valarray representing polynomial coefficients.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::valarray of all roots.
 */
CArry aberth(const DArry &coeffs, const double tol=1e-13);

/**
 * \brief Root-finding function that implements method from Yan & Chieng (2006)
 * 
 * The method proposed by Yan & Chieng can handle roots with multiplicities
 * greater than 1.
 *
 * \param coeffs [in] A std::valarray representing polynomial coefficients.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::valarray of all roots.
 */
CArry yan_and_chieng_2006(const CArry &coeffs, const double tol=1e-13);

/**
 * \brief Root-finding function that implements method from Yan & Chieng (2006)
 * 
 * Overloaded version the accepts std::valarray<double>.
 *
 * \param coeffs [in] A std::valarray representing polynomial coefficients.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::valarray of all roots.
 */
CArry yan_and_chieng_2006(const DArry &coeffs, const double tol=1e-13);

} // end of namespace basic
} // end of namespace simpoly

    
/**
 * \brief Overloaded output stream for CArry.
 *
 * \param os [in] An output stream.
 * \param v [in] A CArry.
 *
 * \return The output stream.
 */
std::ostream &operator<<(std::ostream &os, const simpoly::basic::CArry &v);
    
/**
 * \brief Overloaded output stream for DArry.
 *
 * \param os [in] An output stream.
 * \param v [in] A DArry.
 *
 * \return The output stream.
 */
std::ostream &operator<<(std::ostream &os, const simpoly::basic::DArry &v);
