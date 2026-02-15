add_library(EAM SHARED EAM.cpp)
target_link_libraries(EAM PRIVATE eoncbase ${EON_DEPS})
target_include_directories(EAM PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS EAM)
