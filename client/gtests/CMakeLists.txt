#
# GoogleTests
#
# Created on: 04 Feb 2021
# Author: Rohit Goswami
# Company: University of Iceland
#

# Partially derived from
# https://cliutils.gitlab.io/modern-cmake/chapters/testing/googletest.html

add_subdirectory("${PROJECT_SOURCE_DIR}/subprojects/gprdimer/googletest"
                 "extern/googletest")

mark_as_advanced(
  BUILD_GMOCK
  BUILD_GTEST
  BUILD_SHARED_LIBS
  gmock_build_tests
  gtest_build_samples
  gtest_build_tests
  gtest_disable_pthreads
  gtest_force_shared_crt
  gtest_hide_internal_symbols)

set_target_properties(gtest PROPERTIES FOLDER extern)
set_target_properties(gtest_main PROPERTIES FOLDER extern)
set_target_properties(gmock PROPERTIES FOLDER extern)
set_target_properties(gmock_main PROPERTIES FOLDER extern)

# Test Macros
macro(package_add_test TESTNAME WORKDIR)
  add_executable(${TESTNAME} ${ARGN})
  target_link_libraries(
    ${TESTNAME}
    gtest
    gmock
    gtest_main
    eonlib
    potentials
    gprdimer)
  gtest_discover_tests(
    ${TESTNAME}
    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/gtests/${WORKDIR}"
    PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}/gtests/${WORKDIR}")
  set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
endmacro()

# Colored Verbose Output
add_custom_target(check COMMAND env CTEST_OUTPUT_ON_FAILURE=1 GTEST_COLOR=1
                                ${CMAKE_CTEST_COMMAND})

# Test I
package_add_test(atomconfgpr "data/saddle_search" AtomsConfigurationTest.cpp)
package_add_test(obsgpr "data/saddle_search" ObsTest.cpp)
package_add_test(impdimereon "data/saddle_search" ImpDimerTest.cpp) # pos.con and others are hard-coded
package_add_test(gprdimereon "data/gpr_dimer/morse" GPRDimerTest.cpp) # pos.con and others are hard-coded
