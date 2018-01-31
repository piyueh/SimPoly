/**
 * \file operations.h
 * \brief Operations of polynomials defined by a std::vector<double>.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-29
 */

# pragma once

# include <complex>
# include <valarray>


namespace polynomial
{
namespace op
{
    
    
std::valarray<double> add(
        const std::valarray<double> &p1, const std::valarray<double> &p2);

std::valarray<double> multiply(
        const std::valarray<double> &p1, const std::valarray<double> &p2);

std::valarray<double> divide(
        const std::valarray<double> &p1, const std::valarray<double> &p2);

double evaluate(const double* const &bg, const int n, const double x);

double evaluate(const std::valarray<double> &coeffs, const double x);

std::complex<double> evaluate(const std::complex<double>* const &bg, 
        const int len, const std::complex<double> x);
    
std::complex<double> evaluate(const std::valarray<std::complex<double>> &coeffs, 
        const std::complex<double> x);

double evaluate_root(const double l, 
        const double* const &bg, const int n, const double x);

double evaluate_root(const double l, 
        const std::valarray<double> &roots, const double x);

std::valarray<double> derivative(const std::valarray<double> &coeffs);

std::valarray<double> integral(const std::valarray<double> &coeffs);

std::valarray<double> find_roots(const std::valarray<double> &coeffs, 
        const double tol = 1e-13, 
        const std::valarray<double> &guess=std::valarray<double>(0));

std::valarray<std::complex<double>> find_roots_complex(
        const std::valarray<double> &coeffs, const double tol = 1e-13, 
        const std::valarray<double> &guess=std::valarray<double>(0));

std::valarray<std::complex<double>> find_roots_complex(
        const std::valarray<double> &coeffs, const double tol = 1e-13, 
        const std::valarray<std::complex<double>> &guess=std::valarray<std::complex<double>>(0));

std::valarray<std::complex<double>> find_roots_complex(
        const std::valarray<std::complex<double>> &coeffs, const double tol = 1e-13, 
        const std::valarray<std::complex<double>> &guess=std::valarray<std::complex<double>>(0));

std::valarray<double> find_coefficients(
        const double l, const std::valarray<double> &roots);
    
    
} // end of namespace op
} // end of namespace polynomial
