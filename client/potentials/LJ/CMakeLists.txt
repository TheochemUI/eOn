add_library(LJ SHARED LJ.cpp)
target_link_libraries(LJ PRIVATE eoncbase ${EON_DEPS})
target_include_directories(LJ PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS LJ)
