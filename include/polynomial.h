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
    
class Polynomial
{
    public:
        
        Polynomial() = default;
        Polynomial(Polynomial &&p);
        Polynomial(const basic::DArry &coef);
        Polynomial(const double l, const basic::DArry &rroots);
        Polynomial(const double l, const basic::CArry &croots);
        Polynomial(const double l, 
                const basic::DArry &rroots, const basic::CArry &croots);
        Polynomial(const basic::DArry &coef, const basic::DArry &rroots);
        Polynomial(const basic::DArry &coef, const basic::CArry &croots);
        Polynomial(const basic::DArry &coef, 
                const basic::DArry &roots, const basic::CArry &croots);
        
        virtual ~Polynomial() = default;
        
        void set(const basic::DArry &coef);
        void set(const double l, const basic::DArry &r);
        void set(const double l, const basic::CArry &r);
        void set(const double l, 
                const basic::DArry &rr, const basic::CArry &rc);
        void set(const basic::DArry &coef, const basic::DArry &rroots);
        void set(const basic::DArry &coef, const basic::CArry &croots);
        void set(const basic::DArry &coef, 
                const basic::DArry &roots, const basic::CArry &croots);
        void set(const int d, const double value);
        
        unsigned degree() const;
        basic::DArry coef() const;
        basic::CArry roots() const;
        basic::DArry real_roots() const;
        basic::CArry cmplx_roots() const;
        
        Polynomial deriv() const;
        Polynomial integ() const;
        
        double operator()(const double x) const;
    
    protected:
        
        std::string _type; ///< the type of this polynomial
        char _prefer; ///< preferred to evaluate values; either 'r' or 'c'
        
        unsigned _d; ///< degree of this polynomial
        unsigned _nrr; ///< number of real roots
        unsigned _nrc; ///< number of complex roots
        
        basic::DArry _coef; ///< coefficient array
        basic::DArry _rroots; ///< array holding real roots
        basic::CArry _croots; ///< array holding complex roots
};

} // end of namespace poly
} // end of namespace simpoly
