add_library(VASP SHARED VASP.cpp)
target_link_libraries(VASP PRIVATE eoncbase ${EON_DEPS})
target_include_directories(VASP PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS VASP)
