/**
 * \file basic.h
 * \brief Basic operations of polynomials defined by a std::vector.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-02-06
 */

# pragma once

# include <complex>
# include <vector>


namespace simpoly
{
namespace basic
{

/** \brief Alias of std::vector<...>.
 *
 * \tparam T Basic type of each entry in the std::vector. */
template <typename T> using Arry = std::vector<T>;

/** \brief Alias of std::complex<double>. */
typedef std::complex<double> Cmplx;

/** \brief Alias of std::vector<double>. */
typedef Arry<double> DArry;

/** \brief Alias of std::vector<std::complex<double>>. */
typedef Arry<Cmplx> CArry;

/**
 * \brief Convert a DArry to CArry.
 *
 * \param p [in] A DArry (i.e., std::vector<double>).
 *
 * \return A CArry (i.e., std::vector<std::complex<double>>).
 */
CArry to_CArry(const DArry &p);

/**
 * \brief Get polynomial coefficients by providing roots.
 *
 * \tparam T Base type of entries in the array of roots.
 * \param l [in] The leading coefficient of the polynomial.
 * \param rts [in] An raw pointer to the array of roots.
 * \param len [in] The length of the array of roots.
 *
 * \return A std::vector representing the coefficients.
 */
template <typename T>
Arry<T> to_coefficients(const T &l, const T* const &rts, const int len);

/**
 * \brief Get polynomial coefficients by providing roots.
 *
 * \tparam T Base type of entries in the array of roots.
 * \param l [in] The leading coefficient of the polynomial.
 * \param rts [in] An std::vector of roots.
 *
 * \return A std::vector representing the coefficients.
 */
template <typename T>
Arry<T> to_coefficients(const T &l, const Arry<T> &rts);


/**
 * \brief Addition of two polynomials defined by std::vector.
 *
 * \tparam T Basic type of each entry in the std::vector.
 * \param p1 [in] Coefficients of the first polynomial.
 * \param p2 [in] Coefficients of the second polynomial.
 *
 * \return Resulting polynomial.
 */
template <typename T>
Arry<T> add(const Arry<T> &p1, const Arry<T> &p2);

/**
 * \brief Addition of a polynomials and a constant term.
 *
 * \tparam T Basic type of each entry in the std::vector.
 * \param p1 [in] Coefficients of the polynomial.
 * \param c [in] Constant.
 *
 * \return Resulting polynomial.
 */
template <typename T>
Arry<T> add(const Arry<T> &p, const T &c);

/**
 * \brief Addition of a polynomials and a constant term.
 *
 * \tparam T Basic type of each entry in the std::vector.
 * \param c [in] Constant.
 * \param p1 [in] Coefficients of the polynomial.
 *
 * \return Resulting polynomial.
 */
template <typename T>
Arry<T> add(const T &c, const Arry<T> &p);

/**
 * \brief Substract the 2nd polynomial from the 1st one.
 *
 * \tparam T Basic type of each entry in the std::vector.
 * \param p1 [in] Coefficients of the first polynomial.
 * \param p2 [in] Coefficients of the second polynomial.
 *
 * \return Resulting polynomial.
 */
template <typename T>
Arry<T> substract(const Arry<T> &p1, const Arry<T> &p2);

/**
 * \brief Substract a constant term from a polynomial.
 *
 * \tparam T Basic type of each entry in the std::vector.
 * \param p [in] Coefficients of the polynomial.
 * \param c [in] The constant term.
 *
 * \return Resulting polynomial.
 */
template <typename T>
Arry<T> substract(const Arry<T> &p, const T &c);

/**
 * \brief Substract a polynomial term from a constant.
 *
 * \tparam T Basic type of each entry in the std::vector.
 * \param c [in] The constant term.
 * \param p [in] Coefficients of the polynomial.
 *
 * \return Resulting polynomial.
 */
template <typename T>
Arry<T> substract(const T &c, const Arry<T> &p);

/**
 * \brief Multiplication of two polynomials defined by std::vector.
 *
 * \tparam T Basic type of each entry in the std::vector.
 * \param p1 [in] Coefficients of the first polynomial.
 * \param p2 [in] Coefficients of the second polynomial.
 *
 * \return Resulting polynomial.
 */
template <typename T>
Arry<T> multiply(const Arry<T> &p1, const Arry<T> &p2);

/**
 * \brief Multiplication of a polynomial and a constant.
 *
 * \tparam T Basic type of each entry in the std::vector.
 * \param p [in] Coefficients of the polynomial.
 * \param c [in] The constant term.
 *
 * \return Resulting polynomial.
 */
template <typename T>
Arry<T> multiply(const Arry<T> &p, const T &c);

/**
 * \brief Multiplication of a polynomial and a constant.
 *
 * \tparam T Basic type of each entry in the std::vector.
 * \param c [in] The constant term.
 * \param p [in] Coefficients of the polynomial.
 *
 * \return Resulting polynomial.
 */
template <typename T>
Arry<T> multiply(const T &c, const Arry<T> &p);

/**
 * \brief Division of two polynomials defined by std::vector.
 *
 * \tparam T Basic type of each entry in the std::vector.
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
 * \brief Division of two polynomials defined by std::vector.
 *
 * \tparam T Basic type of each entry in the std::vector.
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
 * \brief Division of a polynomial and a constant.
 *
 * \tparam T Basic type of each entry in the std::vector.
 * \param p [in] Coefficients of the polynomial.
 * \param c [in] The constant.
 *
 * \return Quotient polynomial.
 */
template <typename T>
Arry<T> divide(const Arry<T> &p, const T &c);

/**
 * \brief Find greatest common divisor of two polynomials.
 *
 * \tparam T Base type of the entries in coefficient array.
 * \param p1 [in] A std::vector representing the first polynomial.
 * \param p2 [in] A std::vector representing the second polynomial.
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
 * array of the polynomial is defined by std::vector.
 *
 * \tparam T Basic type of each entry in the std::vector.
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
 * defined by a std::vector of roots and a leading coefficient.
 *
 * \tparam T Basic type of each entry in the root std::vector.
 * \param l [in] The coefficient of the highest-degree term.
 * \param roots [in] A std::vector of roots.
 * \param x [in] The specified location.
 *
 * \return Polynomial value.
 */
template <typename T>
T evaluate_from_root(const T l, const Arry<T> &roots, const T x);


/**
 * \brief Obtain the derived polynomial of a polynomial.
 *
 * \tparam T Basic type of each entry in the root std::vector.
 * \param coeffs [in] A polynomial defined by std::vector.
 *
 * \return Derived polynomial.
 */
template <typename T>
Arry<T> derivative(const Arry<T> &coeffs);

/**
 * \brief Obtain the undeterministic integral of a polynomial.
 *
 * \tparam T Basic type of each entry in the root std::vector.
 * \param coeffs [in] A polynomial defined by std::vector.
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
 * \param coeffs [in] A std::vector representing polynomial coefficients.
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
 * \param coeffs [in] A std::vector representing polynomial coefficients.
 * \param guess [in] A std::vector of initial guess to all roots.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::vector of all roots.
 */
CArry aberth(const CArry &coeffs, const CArry &guess, const double tol=1e-13);

/**
 * \brief Aberth method for root finding.
 *
 * An overloaded version that provided initial guess is a std::vector of
 * real numbers. The output is still a std::vector of complex numbers.
 *
 * \param coeffs [in] A std::vector representing polynomial coefficients.
 * \param guess [in] A std::vector of initial guess to all roots.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::vector of all roots.
 */
CArry aberth(const CArry &coeffs, const DArry &guess, const double tol=1e-13);

/**
 * \brief Aberth method for root finding.
 *
 * An overloaded version that provided coefficients is a std::vector of
 * real numbers. The output is still a std::vector of complex numbers.
 *
 * \param coeffs [in] A std::vector representing polynomial coefficients.
 * \param guess [in] A std::vector of initial guess to all roots.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::vector of all roots.
 */
CArry aberth(const DArry &coeffs, const CArry &guess, const double tol=1e-13);

/**
 * \brief Aberth method for root finding.
 *
 * An overloaded version that both provided coefficients and initial guess are
 * `std::vector`s of real numbers. The output is still a std::vector of
 * complex numbers.
 *
 * \param coeffs [in] A std::vector representing polynomial coefficients.
 * \param guess [in] A std::vector of initial guess to all roots.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::vector of all roots.
 */
CArry aberth(const DArry &coeffs, const DArry &guess, const double tol=1e-13);

/**
 * \brief Aberth method for root finding.
 *
 * An overloaded version of Aberth method that uses default initial guess,
 * instead of user-provided initial guess.
 *
 * \param coeffs [in] A std::vector representing polynomial coefficients.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::vector of all roots.
 */
CArry aberth(const CArry &coeffs, const double tol=1e-13);

/**
 * \brief Aberth method for root finding.
 *
 * An overloaded version of Aberth method that uses default initial guess,
 * instead of user-provided initial guess. Also, the coefficient array is a
 * std::vector of real numbers, not complex numbers.
 *
 * \param coeffs [in] A std::vector representing polynomial coefficients.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::vector of all roots.
 */
CArry aberth(const DArry &coeffs, const double tol=1e-13);

/**
 * \brief Root-finding function that implements method from Yan & Chieng (2006)
 *
 * The method proposed by Yan & Chieng can handle roots with multiplicities
 * greater than 1.
 *
 * \param coeffs [in] A std::vector representing polynomial coefficients.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::vector of all roots.
 */
CArry yan_and_chieng_2006(const CArry &coeffs, const double tol=1e-13);

/**
 * \brief Root-finding function that implements method from Yan & Chieng (2006)
 *
 * Overloaded version the accepts std::vector<double>.
 *
 * \param coeffs [in] A std::vector representing polynomial coefficients.
 * \param tol [in] Tolerance that mimics zero.
 *
 * \return A std::vector of all roots.
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
