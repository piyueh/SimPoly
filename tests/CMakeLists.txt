# =============================================================================
# GOOGLE TESTS
# =============================================================================

# find GTest package
FIND_PACKAGE(GTest REQUIRED)

# add test targets
ADD_EXECUTABLE(operations ${PROJECT_SOURCE_DIR}/tests/test_operations.cpp)
TARGET_INCLUDE_DIRECTORIES(operations PRIVATE ${GTEST_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(operations ${GTEST_BOTH_LIBRARIES} simpoly)

# use GoogleTest
INCLUDE(GoogleTest)
GTEST_DISCOVER_TESTS(operations EXTRA_ARGS --gtest_repeat=10)
