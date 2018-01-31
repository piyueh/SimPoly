# =============================================================================
# libsimpoly
# =============================================================================

# alias to src folder
SET(SRC ${PROJECT_SOURCE_DIR}/src)

# source files
SET(SOURCE
    ${SRC}/operations.cpp)

# SimPoly library
ADD_LIBRARY(simpoly ${SOURCE})

# installation
INSTALL(TARGETS simpoly
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
