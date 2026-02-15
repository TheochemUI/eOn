add_library(ASE SHARED ASE.cpp)
target_link_libraries(ASE PRIVATE eoncbase ${EON_DEPS})
target_include_directories(ASE PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS ASE)
