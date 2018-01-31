# SimPoly -- Simple Polynomial Toolbox

SimPoly is a C++ library for some basic polynomial operations and objects.
Currently it only supports polynomials with real coefficients (though there
are functions that can find complex roots in SimPoly).

SimPoly is not designed for the purpose of high-performance computing. Some
implementations are naive. However, for some high-performance scientific
simulation methods, polynomial operations or calculations are not performance
bottlenecks. Hence a polynomial library with simple usage, like SimPoly, is
useful in this scenario.

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
