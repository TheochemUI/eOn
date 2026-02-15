add_library(
  Water SHARED
  Water.cpp
  ccl.cpp
  potential_base.cpp
  spce_ccl.cpp
  tip4p_ccl.cpp)
target_link_libraries(Water PRIVATE eoncbase ${EON_DEPS})
target_include_directories(Water PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS Water)
