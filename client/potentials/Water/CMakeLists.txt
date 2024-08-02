add_library(
  Water SHARED
  Water.cpp
  ccl.cpp
  potential_base.cpp
  spce_ccl.cpp
  tip4p_ccl.cpp
  )

# Install
install(TARGETS Water)
