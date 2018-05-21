# SimPoly -- Simple Polynomial

SimPoly is a very simple C++ library intended for basic polynomial operations, series, and special polynomials.
Currently, it only supports polynomials with a single variable and real coefficients (though roots can be complex numbers). 

There is already a polynomial module in the Boost library. 
And that module has become a part of C++17 standard. 
The reasons I wrote my own polynomial library:

* I don't like to use Boost just for its polynomial module.
* I don't like the usage of Boost's polynomial module.
* I want something in C++ similar to polynomial objects of Numpy or Scipy.
* Many HPC systems I've used so far have old C++ compilers and are unlikely to have new compilers in the future.

SimPoly is not something computationally efficient. 
Implementations are naive. 
However, for many numerical methods of interest to me, polynomial operations or calculations are not performance bottlenecks. 
So it's fine to me.
SimPoly is not something intended for doing mathematics research.
It's more like a helper for some numerical schemes, such as pseudospectral methods, hp-finite element methods, etc.

Due to the naive implementations, SimPoly may not be numerically stable. 
For ill-conditioned polynomials, it's users' responsibility to take care of rounding errors. 
Again, however, in many scientific computing areas, polynomials are rarely ill-conditioned (I imagine ...). 
So stability issue may not be users' top concern.
At least it's not my concern.

## Features

* `numpy.polynomial.Polynomial`-like usage
    - arithmetic, including division
    - calculus
    - initialize with either coefficients or roots
    - better evaluation if using roots for initialization
* Capability of obtaining accurate roots with multiplicity greater than 1. (Algorithm proposed by Yan & Chieng (2006)[1].)
* Jacobi family polynomials, including Legendre polynomial
* Radau polynomials

## Example code

Usage example of `Polynomial` objects:

```c++
using namespace simpoly;
poly::Polynomial p1({1, 2, 3, 4}); // initialization by coefficients
poly::Polynomial p2(1.0, {1, 2, 3, 4, 5, 6}); // initialization by roots
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

For more operations of `Polynomial` objects, please refer to the header file `include/polynomial.h`

Jacobi polynomials:

```c++
using namespace simpoly;
double alpha=1.0, beta=0.0, degree=4;
poly::Polynomial p = poly::Jacobi(alpha, beta, degree);
```

Legendre polynomials:

```c++
using namespace simpoly;
poly::Polynomial p = poly::Legendre(degree);
```

Left Radau polynomials:

```c++
using namespace simpoly;
poly::Polynomial p = poly::Radau(degree, poly::PolyType::LEFTRADAU);
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

If CMake complains about not finding GTest, please specify the root directory
of GTest through CMake CMD argument: `-DGTEST_ROOT=<path>`.

The following CMake CMD options are supported:

* `CMAKE_BUILD_TYPE`: either `DEBUG`(default) or `RELEASE`
* `BUILD_SHARED_LIBS`: wither `ON` (default) or `OFF`


## Current development

The following features are still in progress (I will implement them when I need them):

* Operation: scaling and shift
* Special polynomial: Lobatto polynomial
* Polynomial series: Lagrange polynomials, Legendre series, Chebychev series, etc
* Quadratures
* Upgrade `double` to `long double`.

Pull requests are welcome.

## Contact

Through GitHub issue section or email pychuang@gwu.edu.

## Reference

[1] Yan, Chang-Dau, and Wei-Hua Chieng. "Method for finding multiple roots of polynomials." Computers & Mathematics with Applications 51, no. 3-4 (2006): 605-620.
