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
enum PolyType: int { General=0, Legendre, Lagrange };

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
     * \param type [in] Type of the polynomial (default: General).
     */
    Polynomial(const basic::DArry &coef, const PolyType &type=General);

    /**
     * \brief Constructor using roots.
     *
     * This will assume all roots of this polynomial lays in real space.
     *
     * \param l [in] The coefficient of the highest degree.
     * \param roots [in] Roots.
     * \param type [in] Type of the polynomial (default: General).
     */
    Polynomial(const double l, const basic::DArry &roots,
            const PolyType &type=General);

    /**
     * \brief Constructor using roots.
     *
     * This will assume all roots of this polynomial lays in complex space.
     *
     * \param l [in] The coefficient of the highest degree.
     * \param roots [in] Roots.
     * \param type [in] Type of the polynomial (default: General).
     */
    Polynomial(const double l, const basic::CArry &roots,
            const PolyType &type=General);

    /**
     * \brief Constructor using roots.
     *
     * We separate real roots and complex roots into two separate vectors.
     *
     * \param l [in] The coefficient of the highest degree.
     * \param rroots [in] Real roots.
     * \param croots [in] Complex roots.
     * \param type [in] Type of the polynomial (default: General).
     */
    Polynomial(const double l, const basic::DArry &rroots,
            const basic::CArry &croots, const PolyType &type=General);

    /**
     * \brief Constructor using coef and roots.
     *
     * This one assumes all roots are in real space. If compiled with RELEASE
     * mode, this function won't check if roots and coefficients match.
     *
     * \param coef [in] Coefficients.
     * \param roots [in] Roots.
     * \param type [in] Type of the polynomial (default: General).
     */
    Polynomial(const basic::DArry &coef, const basic::DArry &roots,
            const PolyType &type=General);

    /**
     * \brief Constructor using coef and roots.
     *
     * This one assumes all roots are in complex space. If compiled with
     * RELEASE mode, this function won't check if roots and coefficients
     * match.
     *
     * \param coef [in] Coefficients.
     * \param roots [in] Roots.
     * \param type [in] Type of the polynomial (default: General).
     */
    Polynomial(const basic::DArry &coef, const basic::CArry &roots,
            const PolyType &type=General);

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
     * \param type [in] Type of the polynomial (default: General).
     */
    Polynomial(const basic::DArry &coef, const basic::DArry &rroots,
            const basic::CArry &croots, const PolyType &type=General);

    /** \brief Destructor. */
    virtual ~Polynomial() = default;

    /**
     * \brief Reset the type of this polynomial.
     *
     * This only resets the type, but coefficients and roots keep unchanged.
     *
     * \param type [in] Type of the polynomial (default: General).
     */
    void set(const PolyType &type);

    /**
     * \brief Reset coefficients and roots using coefficients.
     *
     * \param coef [in] Coefficients.
     * \param tol [in] Tolerance to discard imaginary parts (default: 1e-8).
     */
    void set(const basic::DArry &coef, const double &tol=1e-8);

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
     * \brief Get the number of real roots.
     *
     * \return The number of real roots.
     */
    unsigned n_real_roots() const;

    /**
     * \brief Get the number of complex roots.
     *
     * \return The number of complex roots.
     */
    unsigned n_cmplx_roots() const;

    /**
     * \brief Get a vector of coefficients.
     *
     * \return Coefficients vector.
     */
    basic::DArry coef() const;

    /**
     * \brief Get a vector holding real roots.
     *
     * \return Real roots.
     */
    basic::DArry real_roots() const;

    /**
     * \brief Get a vector holding complex roots.
     *
     * \return complex roots.
     */
    basic::CArry cmplx_roots() const;

    /**
     * \brief Get a complex vector holding all roots, regradless they are real
     *        or complex numbers.
     *
     * \return All roots.
     */
    basic::CArry roots() const;

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
    Polynomial &operator=(const Polynomial &p);

    /**
     * \brief Move assignment.
     *
     * \param p [in] Another Polynomial instance with rvalue type.
     *
     * \return This polynomial with updated content.
     */
    Polynomial &operator=(Polynomial &&p);

protected:

    PolyType _type; ///< the type of this polynomial
    bool _use_roots; ///< indicate if using roots to evaluate values

    unsigned _d; ///< degree of this polynomial
    unsigned _nrr; ///< number of real roots
    unsigned _ncr; ///< number of complex roots

    basic::DArry _coef; ///< coefficient array
    basic::DArry _rroots; ///< array holding real roots
    basic::CArry _croots; ///< array holding complex roots
};

} // end of namespace poly
} // end of namespace simpoly
