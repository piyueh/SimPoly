# =============================================================================
# \file CMakeLists.txt
# \brief CMake settings.
# \author Pi-Yueh Chuang
# \version beta
# \date 2018-01-31
# =============================================================================

CMAKE_MINIMUM_REQUIRED (VERSION 3.10)

# =============================================================================
# DEFAULT VARIABLE VALUES
# =============================================================================

# set default flags for release and debug modes
SET(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")
SET(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined" CACHE STRING "")

# set default build mode
SET(CMAKE_BUILD_TYPE DEBUG CACHE STRING "")


# =============================================================================
# PROJECT SETTINGS
# =============================================================================

# start generator
PROJECT(SimPoly CXX)

# use GNU standard installation directory heirarcy
INCLUDE(GNUInstallDirs)

# check C++11 standard 
SET(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_EXTENSIONS OFF) # use -std=c++11, instead of -std=gnu++11

# other cmake setting
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_COLOR_MAKEFILE ON)

# add options
OPTION(BUILD_SHARED_LIBS "Whether to build shared library." ON)

# project-wide include folders
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)

# SimPoly
ADD_SUBDIRECTORY(src)

# include folder
ADD_SUBDIRECTORY(include)

# enable `make test`
ENABLE_TESTING()
ADD_SUBDIRECTORY(tests)
