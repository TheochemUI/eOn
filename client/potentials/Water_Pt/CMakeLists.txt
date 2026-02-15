add_library(
  Water_Pt SHARED
  Tip4p_Pt.cpp
  zhu_philpott.cpp
  zhu_philpott_parameters.cpp)
target_link_libraries(Water_Pt PRIVATE eoncbase ${EON_DEPS})
target_include_directories(Water_Pt PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS Water_Pt)
