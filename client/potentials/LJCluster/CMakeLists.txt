add_library(LJCluster SHARED LJCluster.cpp)
target_link_libraries(LJCluster PRIVATE eoncbase ${EON_DEPS})
target_include_directories(LJCluster PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS LJCluster)
