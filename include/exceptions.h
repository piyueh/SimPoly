/**
 * \file exceptions.h
 * \brief User-defined exceptions used in Polynomial.
 * \author Pi-Yueh Chuang
 * \version beta
 * \date 2018-01-30
 */

# pragma once

# include <string>
# include <exception>


namespace polynomial
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


class ZeroDegree : public PolynomialErrorGeneral
{
public:
    
    ZeroDegree(const std::string &file, const int &line):
        PolynomialErrorGeneral(file, line, "Zero degree of polynomial detected.") {};
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
                "An array/arrays do not have expected lengths"
                "Length 1: " + std::to_string(l1) +
                "Length 2: " + std::to_string(l2) + ".") {};
    
    const int len1, len2;
};


class InfLoop : public PolynomialErrorGeneral
{
public:
    
    InfLoop(const std::string &file, const int &line):
        PolynomialErrorGeneral(file, line, "Infinite loop.") {};
};


class ComplexRoot : public PolynomialErrorGeneral
{
public:
    
    ComplexRoot(const std::string &file, const int &line,
            const double &r, const double &i):
        real(r), imag(i),
        PolynomialErrorGeneral(file, line, 
                "This polyminal has at least one complex root. Please use "
                "correct function. Root: (" + std::to_string(r) + ", " +
                std::to_string(i) + ").") {};
    
    const double &real, &imag;
};

} // end of namespace exceptions
} // end of namespace polynomial
