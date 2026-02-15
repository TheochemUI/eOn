#
# Catch2 Tests (matching Meson build)
#

# Build the Catch2 amalgamated library + TestMain
add_library(TestMain STATIC
  TestMain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty/catch2/catch_amalgamated.cpp)
target_include_directories(TestMain PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty/catch2
  ${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty
  ${EON_INCLUDE_DIRS}
  ${CMAKE_CURRENT_BINARY_DIR}/..)
target_link_libraries(TestMain PUBLIC eonclib)

# Test macro: creates a Catch2 test executable
macro(eon_add_test TESTNAME WORKDIR)
  add_executable(${TESTNAME}
    ${ARGN}
    ${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty/catch2/catch_amalgamated.cpp)
  target_include_directories(${TESTNAME} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty/catch2
    ${CMAKE_CURRENT_SOURCE_DIR}/../thirdparty
    ${EON_INCLUDE_DIRS}
    ${CMAKE_CURRENT_BINARY_DIR}/..)
  target_link_libraries(${TESTNAME} PRIVATE eonclib)
  if(WORKDIR STREQUAL "")
    add_test(NAME ${TESTNAME} COMMAND ${TESTNAME})
  else()
    add_test(NAME ${TESTNAME} COMMAND ${TESTNAME}
             WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/data/systems/${WORKDIR}")
  endif()
  set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

# Base tests (always built)
eon_add_test(test_zbl "au_si" ZBLPotTest.cpp)
eon_add_test(test_proj_rot_trans "" ProjectOutRotTransTest.cpp)

# Conditional tests
if(WITH_METATOMIC)
  eon_add_test(test_mta "lj38" MetatomicTest.cpp)
endif()

if(WITH_XTB)
  eon_add_test(test_xtb "sulfolene" XTBTest.cpp)
endif()

# Colored verbose output
add_custom_target(check
  COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
  COMMENT "Running tests with verbose output")
