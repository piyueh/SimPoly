# =============================================================================
# libsimpoly
# =============================================================================

# alias to src folder
SET(SRC ${PROJECT_SOURCE_DIR}/src)

# source files
SET(SOURCE
    ${SRC}/basic/operations.cpp
    ${SRC}/basic/evaluations.cpp
    ${SRC}/basic/calculus.cpp
    ${SRC}/basic/root_findings.cpp
    ${SRC}/polynomial/polynomial.cpp
    ${SRC}/polynomial/operators.cpp
    ${SRC}/polynomial/jacobi.cpp
    ${SRC}/polynomial/legendre.cpp
    ${SRC}/polynomial/radau.cpp
    )

# SimPoly library
ADD_LIBRARY(simpoly ${SOURCE})

# installation
INSTALL(TARGETS simpoly
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
