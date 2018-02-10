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

enum PolyType: int { General=0, Legendre, Lagrange };

class Polynomial
{
    public:

        Polynomial() = default;
        Polynomial(const Polynomial &p);
        Polynomial(Polynomial &&p);
        Polynomial(const basic::DArry &coef, const PolyType &type=General);
        Polynomial(const double l, const basic::DArry &roots,
                const PolyType &type=General);
        Polynomial(const double l, const basic::CArry &roots,
                const PolyType &type=General);
        Polynomial(const double l, const basic::DArry &rroots,
                const basic::CArry &croots, const PolyType &type=General);
        Polynomial(const basic::DArry &coef, const basic::DArry &roots,
                const PolyType &type=General);
        Polynomial(const basic::DArry &coef, const basic::CArry &roots,
                const PolyType &type=General);
        Polynomial(const basic::DArry &coef, const basic::DArry &rroots,
                const basic::CArry &croots, const PolyType &type=General);

        virtual ~Polynomial() = default;

        void set(const PolyType &type);
        void set(const basic::DArry &coef, const double &tol=1e-8);
        void set(const double l, const basic::DArry &roots);
        void set(const double l, const basic::CArry &roots);
        void set(const double l,
                const basic::DArry &rroots, const basic::CArry &croots);
        void set(const basic::DArry &coef, const basic::DArry &roots);
        void set(const basic::DArry &coef, const basic::CArry &roots);
        void set(const basic::DArry &coef,
                const basic::DArry &rroots, const basic::CArry &croots);
        void set(const int d, const double value);

        std::string type() const;
        bool use_roots() const;
        unsigned degree() const;
        unsigned n_real_roots() const;
        unsigned c_real_roots() const;
        basic::DArry coef() const;
        basic::DArry real_roots() const;
        basic::CArry cmplx_roots() const;
        basic::CArry roots() const;

        Polynomial deriv() const;
        Polynomial integ() const;

        double operator()(const double x) const;
        Polynomial &operator=(const Polynomial &p);
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
