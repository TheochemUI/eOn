add_library(IMD SHARED IMD.cpp)
target_link_libraries(IMD PRIVATE eoncbase ${EON_DEPS})
target_include_directories(IMD PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS IMD)
