# =============================================================================
# GOOGLE TESTS
# =============================================================================

# find GTest package
FIND_PACKAGE(GTest REQUIRED)

# alias
SET(TEST_SRC ${PROJECT_SOURCE_DIR}/tests)

# add test targets
ADD_EXECUTABLE(operations 
    ${TEST_SRC}/operations/main.cpp
    ${TEST_SRC}/operations/add.cpp
    ${TEST_SRC}/operations/multiply.cpp
    ${TEST_SRC}/operations/divide.cpp
    ${TEST_SRC}/operations/evaluate.cpp)

TARGET_INCLUDE_DIRECTORIES(operations PRIVATE ${GTEST_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(operations ${GTEST_BOTH_LIBRARIES} simpoly)

# use GoogleTest
INCLUDE(GoogleTest)
GTEST_DISCOVER_TESTS(operations EXTRA_ARGS --gtest_repeat=10)
