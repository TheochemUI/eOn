add_library(PyAMFF SHARED PyAMFF.cpp)
target_link_libraries(PyAMFF PRIVATE eoncbase ${EON_DEPS})
target_include_directories(PyAMFF PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS PyAMFF)
