add_library(
  Water_Pt SHARED
  Tip4p_Pt.cpp
  zhu_philpott.cpp
  zhu_philpott_parameters.cpp
  )

# Install
install(TARGETS Water_Pt)
