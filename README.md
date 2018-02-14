# SimPoly -- Simple Polynomial Toolbox

SimPoly is a C++ library for basic polynomial operations, series, and 
special polynomials.
Currently, it only supports polynomials with real coefficients (though roots can be complex numbers). 

There is already a polynomial module in the math module of the Boost library. And that module has become a part of the C++17 standard. To my best knowledge, however, not many compilers have implemented this module in their C++17 implementations. Especially C++ compilers on many HPC systems are still very old. And on the other hand, using Boost library for only to have polynomial operations is overkill. This is one reason I develop my own polynomial library for my future projects. Another reason is that I don't like the usage of the polynomial module in Boost.

SimPoly is not designed for HPC (high-performance computing). Some implementations are naive. However, for many high-performance scientific simulation methods, polynomial operations or calculations are not performance bottlenecks. Hence a polynomial library with simple usage, like SimPoly, is useful in this scenario.

Due to the naive implementations, SimPoly may not be numerically stable. For ill-conditioned polynomials, it's users' responsibility to take care of problems caused by rounding errors. Again, however, in many scientific computing areas, polynomials are rarely ill-conditioned (I imagine ...). So stability issue may not be users' top concern.

## Features
* `numpy.polynomial.Polynomial`-like usage
* Capability of obtaining accurate roots with multiplicity greater than 1. (Algorithm proposed by Yan & Chieng (2006).)

## Example code
```c++
using namespace simpoly;
poly::Polynomial p1({1, 2, 3, 4});
poly::Polynomial p2({1, 2, 3, 4, 5, 6});
poly::Polynomial p3 = p1 * p2;
poly::Polynomial r = p2 % p1;
poly::Polynomial Q, R;
Q = p2.divide(p1, R);
assert(r == R);
p1 += p2;
std::cout << p3 << std::endl;
std::cout << p2.roots() << std::endl;
std::cout << (p2 * p3).interg() << std::endl;
```

## Build and installation

```
$ cd SimPoly
$ mkdir build
$ cmake -DCMAKE_INSTALL_PREFIX=<preferred location> ../
$ make all
$ make test
$ make install
```

The following CMake CMD options are support:

* `CMAKE_BUILD_TYPE`: either `DEBUG`(default) or `RELEASE`
* `BUILD_SHARED_LIBS`: wither `ON` (default) or `OFF`


## Current development
The following features are still in progress:
* Special polynomial: Jacobi-family polynomials, Lagrange base polynomial
* Polynomial series: Lagrange polynomials, Legendre series, Chebychev series, etc
* Quadratures

## Contact
Through GitHub issue section or email pychuang@gwu.edu.
