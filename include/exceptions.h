/**
 * \file exceptions.h
 * \brief User-defined exceptions used in Polynomial.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-30
 */

# pragma once

# include <string>
# include <complex>
# include <exception>

# define __FL__ __FILE__, __LINE__

# ifndef NDEBUG

    # define CHECK_COEFS(c, tol) \
    { \
        if (c.size() == 0) throw simpoly::exceptions::ZeroCoeffsLength(__FL__); \
        if ((c.size() > 1) && (std::abs(c.back()) < tol)) \
            throw simpoly::exceptions::ZeroLeadingCoeff(__FL__); \
    }

# else

    # define CHECK_COEFS(c)

# endif

namespace simpoly
{
namespace exceptions
{

class PolynomialErrorGeneral : public std::exception
{
public:

    PolynomialErrorGeneral(const std::string &file,
            const int &line, const std::string &m):
        msg(m+" (File: "+file+"; Line: "+std::to_string(line)+")") {};

    virtual const char * what() const throw() {return msg.c_str();};

private:

    const std::string msg;
}; // end of class PolynomialErrorGeneral


class DivideByZero : public PolynomialErrorGeneral
{
public:

    DivideByZero(const std::string &file, const int &line):
        PolynomialErrorGeneral(file, line,
                "Divide-by-zero detected.") {};
};


class NegativeDegree : public PolynomialErrorGeneral
{
public:

    NegativeDegree(const std::string &file, const int &line, const int &d):
        degree(d),
        PolynomialErrorGeneral(file, line,
                "Negative degree of polynomial detected. Degree provided is: " +
                std::to_string(d) + ".") {};

    const int degree;
};


class ZeroCoeffsLength : public PolynomialErrorGeneral
{
public:

    ZeroCoeffsLength(const std::string &file, const int &line):
        PolynomialErrorGeneral(file, line,
                "Zero-length coefficient vector/array detected.") {};
};


class ZeroLeadingCoeff : public PolynomialErrorGeneral
{
public:

    ZeroLeadingCoeff(const std::string &file, const int &line):
        PolynomialErrorGeneral(file, line,
                "Zero leading coefficient detected.") {};
};


class NegativeCoeffsLength : public PolynomialErrorGeneral
{
public:

    NegativeCoeffsLength(const std::string &file, const int &line, const int &l):
        len(l),
        PolynomialErrorGeneral(file, line,
                "Negative-length coefficient vector/array detected."
                "Length provided is: " + std::to_string(l) + ".") {};

    const int len;
};


class UnmatchedLength : public PolynomialErrorGeneral
{
public:

    UnmatchedLength(const std::string &file, const int &line,
            const int &l1, const int &l2):
        len1(l1), len2(l2),
        PolynomialErrorGeneral(file, line,
                "An array/arrays do not have expected lengths."
                " Length1: " + std::to_string(l1) +
                " Length2: " + std::to_string(l2) + ".") {};

    const int len1, len2;
};


class InfLoop : public PolynomialErrorGeneral
{
public:

    InfLoop(const std::string &file, const int &line):
        PolynomialErrorGeneral(file, line, "Infinite loop.") {};
};


class FoundComplexNumber : public PolynomialErrorGeneral
{
public:

    FoundComplexNumber(
            const std::string &file, const int &line,
            const std::complex<double> &c):
        number(c),
        PolynomialErrorGeneral(file, line,
                "Found a complex number in a vector that is assumed to be "
                "holding only real numbers. The complex number: (" +
                std::to_string(c.real()) + ", " + std::to_string(c.imag()) +
                ").") {};

    const std::complex<double> &number;
};


class ExpectingZero : public PolynomialErrorGeneral
{
public:

    ExpectingZero(
            const std::string &file, const int &line,
            const double &v):
        value(v),
        PolynomialErrorGeneral(file, line,
                "While expecting a zero, we got a value larger than 1e-12: " +
                std::to_string(v)) {};

    const double &value;
};


class JacobiParameters : public PolynomialErrorGeneral
{
public:

    JacobiParameters(
            const std::string &file, const int &line,
            const double &a, const double &b):
        alpha(a), beta(b),
        PolynomialErrorGeneral(file, line,
            "Alpha or Beta parameter in Jacobi polynomial is wrong: "
            "(Alpha, Beta) = (" +
            std::to_string(a) + ", " + std::to_string(b) + ")") {};

    const double &alpha, &beta;
};

} // end of namespace exceptions
} // end of namespace simpoly
