/**
 * \file polynomial.h
 * \brief Definition of class Polynomial
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-28
 */

# pragma once

# include <memory>

# include "basic.h"

namespace simpoly
{
namespace poly
{

/** \brief Type of polynomial*/
enum PolyType: int { GENERAL=0, JACOBI, LEGENDRE, LAGRANGE };

class Polynomial
{
public:

    /** \brief Default constructor. */
    Polynomial() = default;

    /**
     * \brief Copy constructor.
     *
     * \param p [in] Another Polynomial instance.
     */
    Polynomial(const Polynomial &p);

    /**
     * \brief Move constructor.
     *
     * \param p [in] A rvalue of Polynomial type.
     */
    Polynomial(Polynomial &&p);

    /**
     * \brief Constructor using coefficients.
     *
     * \param coef [in] Coefficients.
     */
    explicit Polynomial(const basic::DArry &coef);

    /**
     * \brief Constructor using roots.
     *
     * This will assume all roots of this polynomial lays in real space.
     *
     * \param l [in] The coefficient of the highest degree.
     * \param roots [in] Roots.
     */
    Polynomial(const double l, const basic::DArry &roots);

    /**
     * \brief Constructor using roots.
     *
     * This will assume all roots of this polynomial lays in complex space.
     *
     * \param l [in] The coefficient of the highest degree.
     * \param roots [in] Roots.
     */
    Polynomial(const double l, const basic::CArry &roots);

    /**
     * \brief Constructor using roots.
     *
     * We separate real roots and complex roots into two separate vectors.
     *
     * \param l [in] The coefficient of the highest degree.
     * \param rroots [in] Real roots.
     * \param croots [in] Complex roots.
     */
    Polynomial(const double l, const basic::DArry &rroots, const basic::CArry &croots);

    /**
     * \brief Constructor using coef and roots.
     *
     * This one assumes all roots are in real space. If compiled with RELEASE
     * mode, this function won't check if roots and coefficients match.
     *
     * \param coef [in] Coefficients.
     * \param roots [in] Roots.
     */
    Polynomial(const basic::DArry &coef, const basic::DArry &roots);

    /**
     * \brief Constructor using coef and roots.
     *
     * This one assumes all roots are in complex space. If compiled with
     * RELEASE mode, this function won't check if roots and coefficients
     * match.
     *
     * \param coef [in] Coefficients.
     * \param roots [in] Roots.
     */
    Polynomial(const basic::DArry &coef, const basic::CArry &roots);

    /**
     * \brief Constructor using coef and roots.
     *
     * This one accepts two separated vecotrs for roots: one for real roots;
     * and the other one for complex roots. If compiled with RELEASE mode,
     * this function won't check if roots and coefficients match.
     *
     * \param coef [in] Coefficients.
     * \param rroots [in] Real roots.
     * \param croots [in] Complex roots.
     */
    Polynomial(const basic::DArry &coef,
            const basic::DArry &rroots, const basic::CArry &croots);

    /** \brief Destructor. */
    virtual ~Polynomial() = default;

    /**
     * \brief Reset the type of this polynomial.
     *
     * This only resets the type, but coefficients and roots keep unchanged.
     *
     * \param type [in] Type of the polynomial (default: GENERAL).
     */
    void set(const PolyType &type);

    /**
     * \brief Reset coefficients and roots using coefficients.
     *
     * \param coef [in] Coefficients.
     */
    void set(const basic::DArry &coef);

    /**
     * \brief Reset coefficients and roots using coefficients.
     *
     * This will assume all roots of this polynomial lays in real space.
     *
     * \param l [in] The coefficient of the highest degree.
     * \param roots [in] Roots.
     */
    void set(const double l, const basic::DArry &roots);

    /**
     * \brief Reset coefficients and roots using coefficients.
     *
     * This will assume all roots of this polynomial lays in complex space.
     *
     * \param l [in] The coefficient of the highest degree.
     * \param roots [in] Roots.
     */
    void set(const double l, const basic::CArry &roots);

    /**
     * \brief Reset coefficients and roots using coefficients.
     *
     * We separate real roots and complex roots into two separate vectors.
     *
     * \param l [in] The coefficient of the highest degree.
     * \param rroots [in] Real roots.
     * \param croots [in] Complex roots.
     */
    void set(const double l, const basic::DArry &rroots, const basic::CArry &croots);

    /**
     * \brief Reset coefficients and roots using coefficients.
     *
     * This one assumes all roots are in real space. If compiled with RELEASE
     * mode, this function won't check if roots and coefficients match.
     *
     * \param coef [in] Coefficients.
     * \param roots [in] Roots.
     */
    void set(const basic::DArry &coef, const basic::DArry &roots);

    /**
     * \brief Reset coefficients and roots using coefficients.
     *
     * This one assumes all roots are in complex space. If compiled with
     * RELEASE mode, this function won't check if roots and coefficients
     * match.
     *
     * \param coef [in] Coefficients.
     * \param roots [in] Roots.
     */
    void set(const basic::DArry &coef, const basic::CArry &roots);

    /**
     * \brief Reset coefficients and roots using coefficients.
     *
     * This one accepts two separated vecotrs for roots: one for real roots;
     * and the other one for complex roots. If compiled with RELEASE mode,
     * this function won't check if roots and coefficients match.
     *
     * \param coef [in] Coefficients.
     * \param rroots [in] Real roots.
     * \param croots [in] Complex roots.
     */
    void set(const basic::DArry &coef,
            const basic::DArry &rroots, const basic::CArry &croots);

    /**
     * \brief Reset coefficients and roots by changing only one coefficient.
     *
     * \param d [in] The target degree to change.
     * \param value [in] New coefficient.
     */
    void set(const int d, const double value);

    /**
     * \brief Get type.
     *
     * \return The type of this polynomial.
     */
    PolyType type() const;

    /**
     * \brief Get a bool indicating whether we are using roots to evaluate
     *        values or not.
     *
     * \return A bool.
     */
    bool use_roots() const;

    /**
     * \brief Get the degree of this polynomial.
     *
     * \return The degree.
     */
    unsigned degree() const;

    /**
     * \brief Get a vector of coefficients.
     *
     * \return Coefficients vector.
     */
    basic::DArry coef() const;

    /**
     * \brief Get the number of real roots.
     *
     * If roots has not been calculated by this instance, then it will claculate
     * roots before returning the number of real roots.
     *
     * \param tol [in] Tolerance mimicing zero used by underlying root-finding
     *        algorithms (default: 1e-12).
     *
     * \return The number of real roots.
     */
    unsigned n_real_roots(const double tol=1e-12) const;

    /**
     * \brief Get the number of complex roots.
     *
     * If roots has not been calculated by this instance, then it will claculate
     * roots before returning the number of complex roots.
     *
     * \param tol [in] Tolerance mimicing zero used by underlying root-finding
     *        algorithms (default: 1e-12).
     *
     * \return The number of complex roots.
     */
    unsigned n_cmplx_roots(const double tol=1e-12) const;

    /**
     * \brief Get a vector holding real roots.
     *
     * If roots has not been calculated by this instance, then it will claculate
     * roots before returning an array of real roots.
     *
     * \param tol [in] Tolerance mimicing zero used by underlying root-finding
     *        algorithms (default: 1e-12).
     *
     * \return Real roots.
     */
    basic::DArry real_roots(const double tol=1e-12) const;

    /**
     * \brief Get a vector holding complex roots.
     *
     * If roots has not been calculated by this instance, then it will claculate
     * roots before returning an array of complex roots.
     *
     * \param tol [in] Tolerance mimicing zero used by underlying root-finding
     *        algorithms (default: 1e-12).
     *
     * \return complex roots.
     */
    basic::CArry cmplx_roots(const double tol=1e-12) const;

    /**
     * \brief Get a complex vector holding all roots, regradless they are real
     *        or complex numbers.
     *
     * If roots has not been calculated by this instance, then it will claculate
     * roots before returning a complex array of all roots.
     *
     * \param tol [in] Tolerance mimicing zero used by underlying root-finding
     *        algorithms (default: 1e-12).
     *
     * \return All roots.
     */
    basic::CArry roots(const double tol=1e-12) const;

    /**
     * \brief Get the derivative of this Polynomial.
     *
     * \return The derivative.
     */
    Polynomial deriv() const;

    /**
     * \brief Get the underministic integral of this Polynomial.
     *
     * \return The integral.
     */
    Polynomial integ() const;

    /**
     * \brief Divide by a polynomial.
     *
     * This is a version of division that returns a quotient and a remainder.
     * For to only get quotient or remainder, please use `quotient()` and
     * `remainder()`.
     *
     * \param divisor [in] Divisor.
     * \param R [out] Remainder.
     *
     * \return Quitient.
     */
    Polynomial divide(const Polynomial &divisor, Polynomial &R);

    /**
     * \brief Divide by a polynomial and return only quotient.
     *
     * Return only the quotient. Any remainder will be discarded. This is equal
     * to operator `/`.
     *
     * \param divisor [in] Divisor.
     *
     * \return Quotient.
     */
    Polynomial quotient(const Polynomial &divisor);

    /**
     * \brief Divide by a polynomial and return only remainder.
     *
     * Return only the quotient. The quotient will be discarded. This is equal
     * to operator `%`.
     *
     * \param divisor [in] Divisor.
     *
     * \return Remainder.
     */
    Polynomial remainder(const Polynomial &divisor);

    /**
     * \brief Overloaded operator().
     *
     * With this operator, to evaluate the value of this polynomial at x, one
     * can use f(x) to get the value.
     *
     * \param x [in] The location to evaluate.
     *
     * \return The value.
     */
    double operator()(const double x) const;

    /**
     * \brief Overloaded operator().
     *
     * With this operator, to evaluate the value of this polynomial at x, one
     * can use f(x) to get the value. This version accepts an array of x.
     *
     * \param x [in] The locations to evaluate.
     *
     * \return The values.
     */
    basic::DArry operator()(const basic::DArry &x) const;

    /**
     * \brief Copy assignment.
     *
     * \param p [in] Another Polynomial instance.
     *
     * \return This polynomial with updated content.
     */
    Polynomial & operator=(const Polynomial &p);

    /**
     * \brief Move assignment.
     *
     * \param p [in] Another Polynomial instance with rvalue type.
     *
     * \return This polynomial with updated content.
     */
    Polynomial & operator=(Polynomial &&p);

    /**
     * \brief Compound operator +=.
     *
     * \param rhs [in] Right hand side Polynomial.
     *
     * \return This polynomial with updated content.
     */
    Polynomial & operator+=(const Polynomial &rhs);

    /**
     * \brief Compound operator +=.
     *
     * \param rhs [in] Right hand side double.
     *
     * \return This polynomial with updated content.
     */
    Polynomial & operator+=(const double &rhs);

    /**
     * \brief Compound operator -=.
     *
     * \param rhs [in] Right hand side Polynomial.
     *
     * \return This polynomial with updated content.
     */
    Polynomial & operator-=(const Polynomial &rhs);

    /**
     * \brief Compound operator -=.
     *
     * \param rhs [in] Right hand side double.
     *
     * \return This polynomial with updated content.
     */
    Polynomial & operator-=(const double &rhs);

    /**
     * \brief Compound operator *=.
     *
     * \param rhs [in] Right hand side Polynomial.
     *
     * \return This polynomial with updated content.
     */
    Polynomial & operator*=(const Polynomial &rhs);

    /**
     * \brief Compound operator *=.
     *
     * \param rhs [in] Right hand side double.
     *
     * \return This polynomial with updated content.
     */
    Polynomial & operator*=(const double &rhs);

    /**
     * \brief Compound operator /=.
     *
     * Only operations with numbers are allowed. For division between
     * polynomials, please use member function `divide`.
     *
     * \param rhs [in] Right hand side double.
     *
     * \return This polynomial with updated content.
     */
    Polynomial & operator/=(const double &rhs);

    /**
     * \brief Comparison operator ==.
     *
     * \param rhs [in] Right hand side Polynomial.
     *
     * \return A bool.
     */
    bool operator==(const Polynomial &rhs) const;

    /**
     * \brief Comparison operator !=.
     *
     * \param rhs [in] Right hand side Polynomial.
     *
     * \return A bool.
     */
    bool operator!=(const Polynomial &rhs) const;

    friend Polynomial operator+(Polynomial lhs, const Polynomial &rhs);
    friend Polynomial operator+(Polynomial lhs, const double &rhs);
    friend Polynomial operator+(const double &lhs, Polynomial rhs);
    friend Polynomial operator-(Polynomial lhs, const Polynomial &rhs);
    friend Polynomial operator-(Polynomial lhs, const double &rhs);
    friend Polynomial operator-(const double &lhs, Polynomial rhs);
    friend Polynomial operator*(Polynomial lhs, const Polynomial &rhs);
    friend Polynomial operator*(Polynomial lhs, const double &rhs);
    friend Polynomial operator*(const double &lhs, Polynomial rhs);
    friend Polynomial operator/(Polynomial lhs, const double &rhs);
    friend Polynomial operator%(Polynomial lhs, const Polynomial &rhs);
    friend std::ostream & operator<<(std::ostream &os, const Polynomial &rhs);
    friend Polynomial divide(const Polynomial &p1,
            const Polynomial &p2, Polynomial &R);
    friend Polynomial quotient(const Polynomial &p1, const Polynomial &p2);
    friend Polynomial remainder(const Polynomial &p1, const Polynomial &p2);

protected:

    PolyType _type; ///< the type of this polynomial
    bool _have_roots; ///< indicate if we already have roots in this instance
    bool _use_roots; ///< indicate if using roots to evaluate values

    unsigned _d; ///< degree of this polynomial
    unsigned _nrr; ///< number of real roots
    unsigned _ncr; ///< number of complex roots

    basic::DArry _coef; ///< coefficient array
    basic::DArry _rroots; ///< array holding real roots
    basic::CArry _croots; ///< array holding complex roots

    /**
     * \brief Underlying private function to get roots.
     *
     * \param tol [in] Tolerance mimics zero (default: 1e-12).
     */
    void _get_roots(const double tol=1e-12) const;
};


/**
 * \brief A factory function creating Jacobi-family polynomials.
 *
 * \param alpha [in] Alpha parameter for Jacobi polynomial.
 * \param beta [in] Beta parameter for Jacobi polynomial.
 * \param n [in] Degree of Jacobi polynomial.
 * \param normalized [in] To make the leading coefficient 1.0 (default: true).
 *
 * \return Jacobi polynomial.
 */
Polynomial Jacobi(const double alpha, const double beta,
        const unsigned n, bool normalized=false);


Polynomial divide(const Polynomial &p1, const Polynomial &p2, Polynomial &R);
Polynomial quotient(const Polynomial &p1, const Polynomial &p2);
Polynomial remainder(const Polynomial &p1, const Polynomial &p2);
Polynomial operator+(Polynomial lhs, const Polynomial &rhs);
Polynomial operator+(Polynomial lhs, const double &rhs);
Polynomial operator+(const double &lhs, Polynomial rhs);
Polynomial operator-(Polynomial lhs, const Polynomial &rhs);
Polynomial operator-(Polynomial lhs, const double &rhs);
Polynomial operator-(const double &lhs, Polynomial rhs);
Polynomial operator*(Polynomial lhs, const Polynomial &rhs);
Polynomial operator*(Polynomial lhs, const double &rhs);
Polynomial operator*(const double &lhs, Polynomial rhs);
Polynomial operator/(Polynomial lhs, const double &rhs);
Polynomial operator%(Polynomial lhs, const Polynomial &rhs);
std::ostream & operator<<(std::ostream &os, const Polynomial &rhs);

} // end of namespace poly
} // end of namespace simpoly
