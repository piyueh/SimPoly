# =============================================================================
# GOOGLE TESTS
# =============================================================================

# find GTest package
FIND_PACKAGE(GTest REQUIRED)

# alias
SET(TEST_SRC ${PROJECT_SOURCE_DIR}/tests)

# add test targets: basic
ADD_EXECUTABLE(basic 
    ${TEST_SRC}/basic/main.cpp
    ${TEST_SRC}/basic/operations1.cpp
    ${TEST_SRC}/basic/operations2.cpp
    ${TEST_SRC}/basic/evaluate.cpp
    ${TEST_SRC}/basic/calculus.cpp
    ${TEST_SRC}/basic/find_roots.cpp
    )

TARGET_INCLUDE_DIRECTORIES(basic PRIVATE ${GTEST_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(basic ${GTEST_BOTH_LIBRARIES} simpoly)

# add test targets: polynomial
ADD_EXECUTABLE(polynomial 
    ${TEST_SRC}/polynomial/main.cpp
    ${TEST_SRC}/polynomial/initialization.cpp
    ${TEST_SRC}/polynomial/reset.cpp
    ${TEST_SRC}/polynomial/evaluation.cpp
    ${TEST_SRC}/polynomial/calculus.cpp
    ${TEST_SRC}/polynomial/operators.cpp
    )

TARGET_INCLUDE_DIRECTORIES(polynomial PRIVATE ${GTEST_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(polynomial ${GTEST_BOTH_LIBRARIES} simpoly)

# use GoogleTest
INCLUDE(GoogleTest)
GTEST_DISCOVER_TESTS(basic EXTRA_ARGS --gtest_repeat=1)
GTEST_DISCOVER_TESTS(polynomial EXTRA_ARGS --gtest_repeat=1)
