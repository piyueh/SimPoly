# =============================================================================
# headers
# =============================================================================

# install headers when using `make install`
INSTALL(DIRECTORY ${PROJECT_SOURCE_DIR}/include
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    FILES_MATCHING PATTERN "*.h")
